{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Rank2Types       #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
module Data.Array.Accelerate.Fourier.Private where

import qualified Data.Array.Accelerate.Convolution.Small as Cyclic
import           Data.Array.Accelerate.Fourier.Sign      (Sign)
import qualified Data.Array.Accelerate.Fourier.Sign      as Sign

import qualified Data.Array.Accelerate.Utility.Sliced    as Sliced
import qualified Data.Array.Accelerate.Utility.Sliced1   as Sliced1

import           Data.Array.Accelerate.LinearAlgebra     (zipExtrudedVectorWith)
import qualified Data.Array.Accelerate.LinearAlgebra     as LinAlg

import           Data.Array.Accelerate.Utility.Lift.Exp  (expr)
import qualified Data.Array.Accelerate.Utility.Lift.Exp  as Exp

import           Data.Array.Accelerate                   ((:.) ((:.)), Acc,
                                                          Array, DIM1, DIM2,
                                                          Elt, Exp, Shape,
                                                          Slice, Z (Z), (!),
                                                          (?))
import qualified Data.Array.Accelerate                   as A
import           Data.Complex                            (Complex ((:+)))


type Transform sh a = Acc (Array sh a) -> Acc (Array sh a)

data SubTransform a =
        SubTransform
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)

data SubTransformPair a =
        SubTransformPair
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)

data SubTransformTriple a =
        SubTransformTriple
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)
           (forall sh. (Shape sh, Slice sh) => Transform (sh:.Int) a)


type PairTransform sh a =
        (Acc (Array sh a), Acc (Array sh a)) ->
        (Acc (Array sh a), Acc (Array sh a))

data SubPairTransform a =
        SubPairTransform
           (forall sh. (Shape sh, Slice sh) => PairTransform (sh:.Int:.Int) a)


cache2 :: (sign ~ Exp (Sign b), a ~ Exp (Complex b), A.RealFloat b, A.Elt (Complex b)) =>
   sign -> a
cache3 :: (sign ~ Exp (Sign b), a ~ Exp (Complex b), A.RealFloat b, A.Elt (Complex b)) =>
   sign -> (a,a)
cache4 :: (sign ~ Exp (Sign b), a ~ Exp (Complex b), A.RealFloat b, A.Elt (Complex b)) =>
   sign -> (a,a,a)
cache5 ::
   (sign ~ Exp (Sign b), a ~ Exp (Complex b),
    A.RealFloat b, A.FromIntegral Int b, A.Elt(Complex b)) =>
   sign -> (a,a,a,a)

cache2 _sign = -1

cache3 sign =
   let sqrt3d2 = sqrt 3 / 2
       mhalf = -1/2
       s = Sign.toSign sign
   in  (A.lift $ mhalf :+ s*sqrt3d2,
        A.lift $ mhalf :+ (-s)*sqrt3d2)

cache4 sign =
   let s = Sign.toSign sign
   in  (A.lift $ 0 :+ s, -1, A.lift $ 0 :+ (-s))

cache5 sign =
   let z = Sign.cisRat sign 5
   in  (z 1, z 2, z 3, z 4)


flatten2 ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array sh (a,a)) ->
   Acc (Array (sh:.Int) a)
flatten2 x =
   A.generate
      (Exp.indexCons (A.shape x) (A.constant 2))
      (Exp.modify (expr :. expr) $
       \(ix :. k)  ->  let xi = x ! ix in k A.== 0 ? (A.fst xi, A.snd xi))

transform2 ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Exp a -> Transform (sh:.Int) a
transform2 z arr =
   flatten2 $
   A.zipWith (\x0 x1 -> A.lift (x0+x1, x0+z*x1))
      (A.slice arr (A.lift $ A.Any :. (0::Int)))
      (A.slice arr (A.lift $ A.Any :. (1::Int)))


flatten3 ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array sh (a,a,a)) ->
   Acc (Array (sh:.Int) a)
flatten3 x =
   A.generate
      (Exp.indexCons (A.shape x) (A.constant (3::Int)))
      (Exp.modify (expr :. expr) $
       \(ix :. k)  ->
          let (x0,x1,x2) = A.unlift $ x ! ix
          in  flip (A.caseof k) x0 $
                 ((A.==1), x1) :
                 ((A.==2), x2) :
                 [])

transform3 ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   (Exp a, Exp a) -> Transform (sh:.Int) a
transform3 (z,z2) arr =
   flatten3 $
   A.zipWith3
      (\x0 x1 x2 ->
         let ((s,_), (zx1,zx2)) = Cyclic.sumAndConvolvePair (x1,x2) (z,z2)
         in  A.lift (x0+s, x0+zx1, x0+zx2))
      (A.slice arr (A.lift $ A.Any :. (0::Int)))
      (A.slice arr (A.lift $ A.Any :. (1::Int)))
      (A.slice arr (A.lift $ A.Any :. (2::Int)))


flatten4 ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array sh (a,a,a,a)) ->
   Acc (Array (sh:.Int) a)
flatten4 x =
   A.generate
      (Exp.indexCons (A.shape x) (A.constant (4::Int)))
      (Exp.modify (expr :. expr) $
       \(ix :. k)  ->
          let (x0,x1,x2,x3) = A.unlift $ x ! ix
          in  flip (A.caseof k) x0 $
                 ((A.==1), x1) :
                 ((A.==2), x2) :
                 ((A.==3), x3) :
                 [])

transform4 ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   (Exp a, Exp a, Exp a) -> Transform (sh:.Int) a
transform4 (z,z2,z3) arr =
   flatten4 $
   A.zipWith4
      (\x0 x1 x2 x3 ->
         let x02a = x0+x2; x02b = x0+z2*x2
             x13a = x1+x3; x13b = x1+z2*x3
         in  A.lift (x02a+   x13a, x02b+z *x13b,
                     x02a+z2*x13a, x02b+z3*x13b))
      (A.slice arr (A.lift $ A.Any :. (0::Int)))
      (A.slice arr (A.lift $ A.Any :. (1::Int)))
      (A.slice arr (A.lift $ A.Any :. (2::Int)))
      (A.slice arr (A.lift $ A.Any :. (3::Int)))


flatten5 ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array sh (a,a,a,a,a)) ->
   Acc (Array (sh:.Int) a)
flatten5 x =
   A.generate
      (Exp.indexCons (A.shape x) (A.constant (5::Int)))
      (Exp.modify (expr :. expr) $
       \(ix :. k)  ->
          let (x0,x1,x2,x3,x4) = A.unlift $ x ! ix
          in  flip (A.caseof k) x0 $
                 ((A.==1), x1) :
                 ((A.==2), x2) :
                 ((A.==3), x3) :
                 ((A.==4), x4) :
                 [])


{-
Use Rader's trick for mapping the transform to a convolution
and apply Karatsuba's trick at two levels (i.e. total three times)
to that convolution.

0 0 0 0 0
0 1 2 3 4
0 2 4 1 3
0 3 1 4 2
0 4 3 2 1

Permutation.T: 0 1 2 4 3

0 0 0 0 0
0 1 2 4 3
0 2 4 3 1
0 4 3 1 2
0 3 1 2 4
-}
transform5 ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   (Exp a, Exp a, Exp a, Exp a) -> Transform (sh:.Int) a
transform5 (z1,z2,z3,z4) arr =
   flatten5 $
   A.zipWith5
      (\x0 x1 x2 x3 x4 ->
         let ((s,_), (d1,d2,d4,d3)) =
                Cyclic.sumAndConvolveQuadruple (x1,x3,x4,x2) (z1,z2,z4,z3)
         in  A.lift (x0+s, x0+d1, x0+d2, x0+d3, x0+d4))
      (A.slice arr (A.lift $ A.Any :. (0::Int)))
      (A.slice arr (A.lift $ A.Any :. (1::Int)))
      (A.slice arr (A.lift $ A.Any :. (2::Int)))
      (A.slice arr (A.lift $ A.Any :. (3::Int)))
      (A.slice arr (A.lift $ A.Any :. (4::Int)))


twist ::
   (Shape sh, Slice sh, Elt a) =>
   Exp Int -> Transform (sh:.Int:.Int) a
twist fac x =
   let sh :. m :. n = Exp.unlift (expr :. expr :. expr) $ A.shape x
   in  A.backpermute
          (A.lift $ sh :. fac*m :. div n fac)
          (Exp.modify (expr :. expr :. expr) $
           \(globalIx :. k :. j) -> globalIx :. div k fac :. fac*j + mod k fac)
          x


merge ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array (sh:.Int:.Int) a) ->
   Acc (Array (sh:.Int) a)
merge x =
   let sh :. m :. n = Exp.unlift (expr :. expr :. expr) $ A.shape x
   in  A.backpermute
          (A.lift $ sh :. m*n)
          (Exp.modify (expr :. expr) $
           \(ix :. k) -> ix :. mod k m :. div k m)
          x

stack ::
   (Shape sh, Slice sh, Elt a) =>
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int:.Int) a)
stack x y =
   A.generate
      (Exp.modify (expr :. expr)
         (\(sh :. n) -> sh :. (2::Int) :. n)
         (A.shape x))
      (Exp.modify (expr :. expr :. expr) $
       \(globalIx :. evenOdd :. k) ->
          let ix = A.lift $ globalIx :. k
          in  evenOdd A.== 0 ? (x ! ix, y ! ix))


{- |
twiddle factors for radix-2 Cooley-Tukey transforms
-}
twiddleFactors2 ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Acc (A.Vector (Complex a))
twiddleFactors2 sign len2 =
   A.generate (A.lift $ Z:.len2) $ twiddle2 sign len2 . A.indexHead

twiddle2 ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Exp (Complex a)
twiddle2 sign n2i ki =
   let n2 = A.fromIntegral n2i
       k  = A.fromIntegral ki
   in  Sign.cis sign $ pi*k/n2


twiddleFactors ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Acc (Array DIM2 (Complex a))
twiddleFactors sign lenk lenj =
   A.generate (A.lift $ Z:.lenk:.lenj) $
   Exp.modify (expr :. expr :. expr) $
   \(_z :. k :. j) -> twiddle sign (lenk*lenj) k j

twiddle ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Exp Int -> Exp (Complex a)
twiddle sign n k j =
   Sign.cisRat sign n $ mod (k*j) n


transformRadix2InterleavedTime ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Acc (Array DIM1 a) ->
   Transform (sh:.Int:.Int) a ->
   Transform (sh:.Int) a
transformRadix2InterleavedTime twiddles subTransform arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
       len2 = div len 2
       subs =
          subTransform $
          if True
            then Sliced.sliceHorizontal (A.lift $ Z:.(2::Int):.len2) arr
            else
               LinAlg.transpose $
               A.reshape (A.lift $ sh:.len2:.(2::Int)) arr
       evens = A.slice subs (A.lift $ A.Any :. (0::Int) :. A.All)
       odds =
          zipExtrudedVectorWith (*) twiddles $
          A.slice subs (A.lift $ A.Any :. (1::Int) :. A.All)
   in  A.zipWith (+) evens odds  A.++  A.zipWith (-) evens odds


initSplitRadix ::
   (Slice sh, Shape sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Acc (Array (sh:.Int) a) ->
   (Acc (Array (sh:.Int:.Int) a), Acc (Array (sh:.Int:.Int) a))
initSplitRadix arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
   in  (A.replicate (A.lift $ A.Any :. (1::Int) :. A.All) arr,
        A.fill (A.lift $ sh:.(0::Int):.div len 2) 0)

finishSplitRadix ::
   (Slice sh, Shape sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Acc (Array (sh:.Int:.Int) a) -> Acc (Array (sh:.Int) a)
finishSplitRadix =
   flip A.slice (A.lift $ A.Any :. (0::Int) :. A.All)


initSplitRadixFlat ::
   (Slice sh, Shape sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Acc (Array (sh:.Int) a) ->
   (Acc (Array DIM2 a), Acc (Array DIM2 a))
initSplitRadixFlat arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
   in  (A.reshape (A.lift $ Z :. A.shapeSize sh :. len) arr,
        A.fill (A.lift $ Z:.(0::Int):.div len 2) 0)

finishSplitRadixFlat ::
   (Slice sh, Shape sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Exp (sh:.Int) -> Acc (Array DIM2 a) -> Acc (Array (sh:.Int) a)
finishSplitRadixFlat = A.reshape


imagSplitRadixPlain ::
   (Num a) =>
   Sign a -> Complex a
imagSplitRadixPlain sign = 0 :+ Sign.getSign sign

imagSplitRadix ::
   (A.Num a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp (Complex a)
imagSplitRadix sign =
   A.lift (0 :+ Sign.toSign sign)

ditSplitRadixReorder ::
   (Slice sh, Shape sh, Elt a) =>
   PairTransform (sh:.Int:.Int) a
ditSplitRadixReorder (arr2, arr1) =
   let evens = Sliced.sieve 2 0 arr2
       odds  = Sliced.sieve 2 1 arr2
   in  (Sliced1.append evens arr1, twist 2 odds)

ditSplitRadixBase ::
   (Slice sh, Shape sh, A.RealFloat a, A.Elt (Complex a)) =>
   PairTransform (sh:.Int:.Int) (Complex a)
ditSplitRadixBase (arr2, arr1) = (transform2 (-1) arr2, arr1)

ditSplitRadixStep ::
   (Slice sh, Shape sh, a ~ Complex b, A.RealFloat b, A.Elt (Complex b)) =>
   Exp a ->
   (Acc (Array DIM1 a), Acc (Array DIM1 a)) ->
   PairTransform (sh:.Int:.Int) a
ditSplitRadixStep imag (twiddles1, twiddles3) (u, zIntl) =
   let twiddledZEven =
          zipExtrudedVectorWith (*) twiddles1 $ Sliced1.sieve 2 0 zIntl
       twiddledZOdd =
          zipExtrudedVectorWith (*) twiddles3 $ Sliced1.sieve 2 1 zIntl
       zSum = A.zipWith (+) twiddledZEven twiddledZOdd
       zDiff = A.map (imag *) $ A.zipWith (-) twiddledZEven twiddledZOdd
       zComplete = zSum A.++ zDiff
   in  (A.zipWith (+) u zComplete
        A.++
        A.zipWith (-) u zComplete,
           Sliced1.drop (Sliced1.length zComplete) u)


twiddleSR ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Exp Int -> Exp (Complex a)
twiddleSR sign n4i ki ji =
   let n4 = A.fromIntegral n4i
       k  = A.fromIntegral ki
       j  = A.fromIntegral ji
   in  Sign.cis sign $ pi*(k*j)/(2*n4)

twiddleFactorsSR ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Acc (Array DIM1 (Complex a))
twiddleFactorsSR sign len4 k =
   A.generate (A.lift $ Z:.len4) $ twiddleSR sign len4 k . A.indexHead

twiddleFactorsSRPair ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int ->
   (Acc (Array DIM1 (Complex a)), Acc (Array DIM1 (Complex a)))
twiddleFactorsSRPair sign len4 =
   (twiddleFactorsSR sign len4 1,
    twiddleFactorsSR sign len4 3)


cycleDim3 ::
   (Slice sh, Shape sh, Elt a) =>
   Acc (Array (sh:.Int:.Int:.Int) a) ->
   Acc (Array (sh:.Int:.Int:.Int) a)
cycleDim3 arr =
   A.backpermute
      (Exp.modify (expr:.expr:.expr:.expr)
         (\(sh:.k:.m:.n) -> (sh:.n:.k:.m)) $
       A.shape arr)
      (Exp.modify (expr:.expr:.expr:.expr)
         (\(ix:.n:.k:.m) -> (ix:.k:.m:.n)))
      arr


chirp ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp a -> A.Acc (A.Array DIM1 (Complex a))
chirp sign padLen lenFloat =
   A.generate (A.index1 padLen) $
   \ix ->
      let k = A.unindex1 ix
          sk = A.fromIntegral (padLen A.> 2*k ? (k, k-padLen))
      in  Sign.cis sign (pi*sk*sk/lenFloat)
