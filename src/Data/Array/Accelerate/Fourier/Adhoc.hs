{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
{- |
The implementations in this module work entirely in the 'A.Acc' domain.
This means that they can be applied to any array
without knowing their extent on the Haskell side.
The downside is, that they cannot share any preprocessing.
-}
module Data.Array.Accelerate.Fourier.Adhoc (
   Transform,
   transform,

   ditSplitRadix,
   dit2,
   dit235,

   ceiling5Smooth,
   transformChirp2,
   transformChirp235,

   Sign,
   forward,
   inverse,

   transform2d,
   transform3d,

   SubTransform(SubTransform),
   ) where

import           Data.Array.Accelerate.Fourier.Private  (SubTransform (SubTransform),
                                                         SubTransformPair (SubTransformPair),
                                                         Transform, twist)
import qualified Data.Array.Accelerate.Fourier.Private  as Fourier
import           Data.Array.Accelerate.Fourier.Sign     (Sign)
import qualified Data.Array.Accelerate.Fourier.Sign     as Sign
import           Data.Array.Accelerate.Fourier.Utility  (scaleDown)

import           Data.Array.Accelerate.LinearAlgebra    (zipExtrudedMatrixWith,
                                                         zipExtrudedVectorWith)
import qualified Data.Array.Accelerate.LinearAlgebra    as LinAlg

import qualified Data.Array.Accelerate.Utility.Arrange  as Arrange
import           Data.Array.Accelerate.Utility.Lift.Acc (acc)
import qualified Data.Array.Accelerate.Utility.Lift.Acc as Acc
import           Data.Array.Accelerate.Utility.Lift.Exp (expr)
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import           Data.Array.Accelerate.Utility.Ord      (argminimum)
import qualified Data.Array.Accelerate.Utility.Sliced   as Sliced
import qualified Data.Array.Accelerate.Utility.Sliced1  as Sliced1

import           Data.Array.Accelerate.Data.Complex     (Complex)
import qualified Data.Array.Accelerate.Data.Complex     as Complex

import           Data.Array.Accelerate                  ((:.) ((:.)), DIM1, Exp,
                                                         Shape, Slice, Z (Z))
import qualified Data.Array.Accelerate                  as A
import           Data.Array.Accelerate.Data.Bits        ((.&.))



forward, inverse :: (A.Num a) => Exp (Sign a)
forward = Sign.forwardExp
inverse = Sign.inverseExp


{- |
Automatically choose transformation algorithms
according to the size of the array.
However, they are not as sophisticated as the algorithms in
"Data.Array.Accelerate.Fourier.Planned".
-}
transform ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int) (Complex a)
transform sign arr =
   let len = Sliced.length arr
   in  A.acond (len A.<= 1) arr $
          let (pow2, smooth5) = is2or5smooth len
          in  A.acond pow2 (ditSplitRadixLoop sign arr) $
              A.acond smooth5 (dit235 sign arr) $
              transformChirp235 sign arr

is2or5smooth :: Exp Int -> (Exp Bool, Exp Bool)
is2or5smooth len =
   let maxPowerOfTwo = len .&. negate len
       lenOdd = div len maxPowerOfTwo
   in  (lenOdd A.== 1,
        (divideMaxPower 5 $ divideMaxPower 3 lenOdd) A.== 1)

divideMaxPower :: Exp Int -> Exp Int -> Exp Int
divideMaxPower fac =
   A.while (\n -> mod n fac A.== 0) (flip div fac)


{- |
Split-Radix for power-of-two sizes.
-}
ditSplitRadix ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int) (Complex a)
ditSplitRadix sign arr =
   A.acond
      (Sliced.length arr A.<= 1)
      arr (ditSplitRadixLoop sign arr)

ditSplitRadixLoop ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int) (Complex a)
ditSplitRadixLoop sign =
   Fourier.finishSplitRadix . A.afst
   .
   A.awhile
      (\x -> A.unit $ (Sliced1.length $ A.asnd x) A.> 0)
      (Acc.modify (acc, acc) $
       \(arr2, arr1) ->
         Fourier.ditSplitRadixStep
            (Fourier.imagSplitRadix sign)
            (Fourier.twiddleFactorsSRPair sign $ Sliced.length arr1)
            (arr2, arr1))
   .
   Acc.modify (acc, acc) Fourier.ditSplitRadixBase
   .
   A.awhile
      (\x -> A.unit $ (Sliced.length $ A.asnd x) A.> 1)
      (Acc.modify (acc, acc) Fourier.ditSplitRadixReorder)
   .
   A.lift . Fourier.initSplitRadix


{- |
Decimation in time for power-of-two sizes.
-}
dit2 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int) (Complex a)
dit2 sign =
   flip A.slice (A.lift $ A.Any :. (0::Int) :. A.All)
   .
   A.awhile
      (\x -> A.unit $ Sliced1.length x A.> 1)
      (ditStep sign)
   .
   A.awhile
      (\x -> A.unit $ Sliced.length x A.> 1)
      (twist 2)
   .
   A.replicate (A.lift $ A.Any :. (1::Int) :. A.All)

ditStep ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int:.Int) (Complex a)
ditStep sign x =
   let twiddles = Fourier.twiddleFactors2 sign $ Sliced.length x
       evens = Sliced1.sieve 2 0 x
       odds = zipExtrudedVectorWith (*) twiddles $ Sliced1.sieve 2 1 x
   in  A.zipWith (+) evens odds  A.++  A.zipWith (-) evens odds


{- |
Decimation in time for sizes that are composites of the factors 2, 3 and 5.
These sizes are known as 5-smooth numbers or the Hamming sequence.
<http://oeis.org/A051037>.
-}
dit235 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Transform (sh:.Int) (Complex a)
dit235 sign =
   flip A.slice (A.lift $ A.Any :. (0::Int) :. A.All)
   .
   A.afst
   .
   A.awhile
      (\x -> A.unit $ (A.length $ A.asnd x) A.> 0)
      (Acc.modify (acc,acc) $
       \(arr,factors) ->
         let fac = factors A.! A.index1 0
         in  (dit235Step sign fac arr, Sliced.tail factors))
   .
   A.awhile
      (\x -> A.unit $ (Sliced.length $ A.afst x) A.> 1)
      (Acc.modify (acc,acc) $
       \(arr,factors) ->
         let divides k n = mod n k A.== 0
             caseFactor k = (divides k, k)
             len = Sliced.length arr
             factor =
                flip (A.caseof len) 2 $
                   caseFactor 3 :
                   caseFactor 4 :
                   caseFactor 5 :
                   []
         in  (twist factor arr, Sliced.consExp factor factors))
   .
   A.lift . flip (,) (A.fill (A.index1 0) 0)
   .
   A.replicate (A.lift $ A.Any :. (1::Int) :. A.All)

dit235Step ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Exp Int ->
   Transform (sh:.Int:.Int) (Complex a)
dit235Step sign fac x =
   let (sh:.count:.len) = Exp.unlift (expr:.expr:.expr) $ A.shape x
       twiddled =
          LinAlg.transpose .
          zipExtrudedMatrixWith (*) (Fourier.twiddleFactors sign fac len) .
          A.reshape (A.lift (sh :. div count fac :. fac :. len))
          $
          x
   in  Fourier.merge $
       A.acond (fac A.== 5) (Fourier.transform5 (Fourier.cache5 sign) twiddled) $
       A.acond (fac A.== 4) (Fourier.transform4 (Fourier.cache4 sign) twiddled) $
       A.acond (fac A.== 3) (Fourier.transform3 (Fourier.cache3 sign) twiddled) $
       Fourier.transform2 (Fourier.cache2 sign) twiddled


{- |
Next greater or equal 5-smooth number as needed by 'dit235'.
-}
ceiling5Smooth :: Exp Int -> Exp Int
ceiling5Smooth n =
   Exp.modify (expr,expr,expr)
      (\(e2, e3, e5) -> pow e2 2 * pow e3 3 * pow e5 5) $
   A.snd $ ceiling5SmoothFloat $
   (A.fromIntegral n :: Exp Double)

ceiling5SmoothFloat ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Exp a -> Exp (a, (Int, Int, Int))
ceiling5SmoothFloat n =
   let d3 = ceilingLogBase 3 n
       d5 = ceilingLogBase 5 n
   in  A.the $ argminimum $
       A.generate (A.lift $ Z:.d5:.d3) $
       Exp.modify (expr:.expr:.expr) $
       \(_z:.e5:.e3) ->
          let p53 = 5 ** A.fromIntegral e5 * 3 ** A.fromIntegral e3
              e2 = max 0 $ ceilingLogBase 2 $ n/p53
          in  (p53 * 2 ** A.fromIntegral e2, (e2, e3, e5))

{-
Should be more efficient than ceiling5SmoothFloat,
but sometimes misses optimal results due to rounding errors.
-}
_ceiling5SmoothLog ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Exp a -> Exp (a, (Int, Int, Int))
_ceiling5SmoothLog n =
   let log3 = logBase 2 3
       log5 = logBase 2 5
       logN = logBase 2 n
       d3 = max 1 $ A.ceiling $ logN / log3
       d5 = max 1 $ A.ceiling $ logN / log5
   in  A.the $ argminimum $
       A.generate (A.lift $ Z:.d5:.d3) $
       Exp.modify (expr:.expr:.expr) $
       \(_z:.e5:.e3) ->
          let logP53 = log5 * A.fromIntegral e5 + log3 * A.fromIntegral e3
              e2 = max 0 $ A.ceiling $ logN-logP53
          in  (logP53 + A.fromIntegral e2, (e2, e3, e5))

_ceiling5SmoothFloat ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Exp a -> Exp (a, (Int, Int, Int))
_ceiling5SmoothFloat n =
   let powers base =
          A.scanl (*) 1 $
          A.fill (A.index1 $ ceilingLogBase base n) $
          A.fromIntegral base
   in  A.the $ argminimum $
       Arrange.mapWithIndex
          (Exp.modify2 (expr:.expr:.expr) expr $
           \(_z:.e5:.e3) p53 ->
              let e2 = max 0 $ ceilingLogBase 2 $ n/p53
              in  (p53 * 2 ** A.fromIntegral e2, (e2, e3, e5))) $
       LinAlg.outer (powers 5) (powers 3)

ceilingLogBase ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Exp Int -> Exp a -> Exp Int
ceilingLogBase base x =
   A.ceiling $ logBase (A.fromIntegral base) x

pow :: Exp Int -> Exp Int -> Exp Int
pow e n = A.the $ A.product $ A.fill (A.index1 e) n


_transformChirp ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Exp Int ->
   (Transform DIM1 (Complex a),
    Transform (sh:.Int) (Complex a),
    Transform (sh:.Int) (Complex a)) ->
   Transform (sh:.Int) (Complex a)
_transformChirp sign padLen (analysis1,analysis,synthesis) arr =
   let len = Sliced.length arr
       chirp = Fourier.chirp sign padLen $ A.fromIntegral len
   in  A.acond (len A.<= 1) arr $
       Sliced.take len $ scaleDown $
       LinAlg.zipExtrudedVectorWith (*) chirp $ synthesis $
       LinAlg.zipExtrudedVectorWith (*)
          (analysis1 $ A.map Complex.conjugate chirp)
          (analysis $ Sliced.pad 0 padLen $
           LinAlg.zipExtrudedVectorWith (*) chirp arr)

transformChirp ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Exp Int ->
   SubTransformPair (Complex a) ->
   Transform (sh:.Int) (Complex a)
transformChirp sign padLen (SubTransformPair analysis synthesis) arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
       chirp = Fourier.chirp sign padLen $ A.fromIntegral len
       spectrum =
          analysis $
          Sliced1.cons
             (A.map Complex.conjugate chirp)
             (A.reshape (A.lift $ A.index1 (A.shapeSize sh) :. padLen) $
              Sliced.pad 0 padLen $
              LinAlg.zipExtrudedVectorWith (*) chirp arr)
   in  A.acond (len A.<= 1) arr $
       Sliced.take len $ scaleDown $
       LinAlg.zipExtrudedVectorWith (*) chirp $ synthesis $
       LinAlg.zipExtrudedVectorWith (*)
          (Sliced1.head spectrum)
          (A.reshape (A.lift $ sh:.padLen) $ Sliced1.tail spectrum)

{- |
Transformation of arbitrary length based on Bluestein on a power-of-two size.
-}
transformChirp2 ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Transform (sh:.Int) (Complex a)
transformChirp2 sign arr =
   transformChirp sign
      (let n = Sliced.length arr
       in  pow (ceilingLogBase 2 (A.fromIntegral (2*n-1) :: Exp Double)) 2)
      (SubTransformPair (ditSplitRadix forward) (ditSplitRadix inverse))
      arr

{- |
Transformation of arbitrary length based on Bluestein on a 5-smooth size.
-}
transformChirp235 ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Transform (sh:.Int) (Complex a)
transformChirp235 sign arr =
   transformChirp sign
      (ceiling5Smooth (2 * Sliced.length arr))
      (SubTransformPair (dit235 forward) (dit235 inverse))
      arr


transform2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   SubTransform (Complex a) ->
   Transform (sh:.Int:.Int) (Complex a)
transform2d (SubTransform trans) =
   LinAlg.transpose . trans .
   LinAlg.transpose . trans

transform3d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   SubTransform (Complex a) ->
   Transform (sh:.Int:.Int:.Int) (Complex a)
transform3d (SubTransform trans) =
   Fourier.cycleDim3 . trans .
   Fourier.cycleDim3 . trans .
   Fourier.cycleDim3 . trans
