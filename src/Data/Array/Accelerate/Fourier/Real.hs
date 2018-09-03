{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
{- |
Compute transforms on real data based on complex-valued transforms.
-}
module Data.Array.Accelerate.Fourier.Real (
   toSpectrum,
   fromSpectrum,

   twoToSpectrum,
   twoToSpectrum2d,
   untangleSpectra,
   untangleSpectra2d,
   untangleCoefficient,

   twoFromSpectrum,
   twoFromSpectrum2d,
   entangleSpectra,
   entangleSpectra2d,
   entangleCoefficient,
   ) where

import qualified Data.Array.Accelerate.Cyclic           as Cyclic
import qualified Data.Array.Accelerate.Fourier.Private  as Fourier
import qualified Data.Array.Accelerate.Fourier.Sign     as Sign

import           Data.Array.Accelerate.Utility.Lift.Exp (expr)
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import qualified Data.Array.Accelerate.Utility.Sliced   as Sliced

import           Data.Array.Accelerate.LinearAlgebra    (zipExtrudedVectorWith)

import           Data.Array.Accelerate.Data.Complex     (Complex ((:+)))
import qualified Data.Array.Accelerate.Data.Complex     as Complex

import           Data.Array.Accelerate                  ((:.) ((:.)), Acc,
                                                         Array, Exp, Shape,
                                                         Slice, (!), (?))
import qualified Data.Array.Accelerate                  as A


{- |
Perform a real-to-complex transform
using a complex-to-complex transform of half size.
Input must have an even size.
Result has the same size as the input, i.e. it is not halved.
-}
toSpectrum ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int) (Complex a) ->
   Acc (Array (sh:.Int) a) -> Acc (Array (sh:.Int) (Complex a))
toSpectrum subTrans arr =
   let n2 = div (Sliced.length arr) 2
       x = subTrans $ complexDeinterleave arr
       xp = A.map (/2) $ A.zipWith (+) x (Cyclic.reverse x)
       xm = A.map (/2) $ A.zipWith (-) x (Cyclic.reverse x)
       twiddles =
          A.map (imagUnit*) $
          Fourier.twiddleFactors2 Sign.forwardExp n2
       evens =
          A.zipWith
             (\xpk xmk -> A.lift $ Complex.real xpk :+ Complex.imag xmk) xp xm
       odds  =
          zipExtrudedVectorWith (*) twiddles $
          A.zipWith
             (\xpk xmk -> A.lift $ Complex.real xmk :+ Complex.imag xpk) xp xm
   in  A.zipWith (-) evens odds
       A.++
       A.zipWith (+) evens odds

complexDeinterleave ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int) a) -> Acc (Array (sh:.Int) (Complex a))
complexDeinterleave arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
   in  A.generate
          (A.lift $ sh :. div len 2)
          (Exp.modify (expr:.expr) $
           \(ix:.j) ->
             arr ! A.lift (ix:.2*j)
             :+
             arr ! A.lift (ix:.2*j+1))


{- |
Perform a complex-to-real transform
using a complex-to-complex of half size.
Input must be self-adjoint and must have an even size.
Result has the same size as the input, i.e. it is not doubled.
-}
fromSpectrum ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int) (Complex a) ->
   Acc (Array (sh:.Int) (Complex a)) -> Acc (Array (sh:.Int) a)
fromSpectrum subTrans spec =
   let n2 = div (Sliced.length spec) 2
       twiddles =
          A.map (imagUnit*) $
          Fourier.twiddleFactors2 Sign.inverseExp n2
       part0 = Sliced.take n2 spec
       part1 = Sliced.drop n2 spec
       fe = A.zipWith (+) part0 part1
       fo =
          zipExtrudedVectorWith (*) twiddles $
          A.zipWith (-) part0 part1
   in  complexInterleave $ subTrans $ A.zipWith (+) fe fo

complexInterleave ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int) (Complex a)) -> Acc (Array (sh:.Int) a)
complexInterleave arr =
   let (sh:.len) = Exp.unlift (expr:.expr) $ A.shape arr
   in  A.generate
          (A.lift $ sh :. 2*len)
          (Exp.modify (expr:.expr) $
           \(ix:.j) ->
             let k = div j 2
                 r = mod j 2
                 x = arr ! A.lift (ix:.k)
             in  r A.== 0 ? (Complex.real x, Complex.imag x))


{- |
Perform a real-to-complex transform of two real inputs
using a complex-to-complex transform of the same size.
Input can have arbitrary size.
-}
twoToSpectrum ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int) (Complex a) ->
   Acc (Array (sh:.Int) (a,a)) ->
   Acc (Array (sh:.Int) (Complex a, Complex a))
twoToSpectrum subTrans =
   untangleSpectra . subTrans .
   A.map (Exp.modify (expr,expr) $ uncurry (:+))

twoToSpectrum2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int:.Int) (Complex a) ->
   Acc (Array (sh:.Int:.Int) (a,a)) ->
   Acc (Array (sh:.Int:.Int) (Complex a, Complex a))
twoToSpectrum2d subTrans =
   untangleSpectra2d . subTrans .
   A.map (Exp.modify (expr,expr) $ uncurry (:+))

{- |
You can transform two real data sets using one complex transform.
This function can be used to untangle the resulting spectrum.
-}
{-
Let f and g be two real valued images.
The spectrum of f+i*g is spec f + i * spec g.
Let 'flip' be the spectrum with negated indices modulo image size.
It holds: flip (spec f) = conj (spec f).

(a + conj b) / 2
  = (spec (f+i*g) + conj (flip (spec (f+i*g)))) / 2
  = (spec f + i*spec g + conj (flip (spec f)) + conj (flip (spec (i*g)))) / 2
  = (2*spec f + i*spec g + conj (i*flip (spec g))) / 2
  = (2*spec f + i*spec g - i * conj (flip (spec g))) / 2
  = spec f

(a - conj b) * (-i/2)
  = (-i*a + conj (-i*b)) / 2
  -> this swaps role of f and g in the proof above
-}
untangleSpectra ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int) (Complex a)) ->
   Acc (Array (sh:.Int) (Complex a, Complex a))
untangleSpectra spec =
   A.zipWith untangleCoefficient spec (Cyclic.reverse spec)

untangleSpectra2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int:.Int) (Complex a)) ->
   Acc (Array (sh:.Int:.Int) (Complex a, Complex a))
untangleSpectra2d spec =
   A.zipWith untangleCoefficient spec (Cyclic.reverse2d spec)

untangleCoefficient ::
   (A.RealFloat a, A.Elt (Complex a)) =>
   Exp (Complex a) -> Exp (Complex a) -> Exp (Complex a, Complex a)
untangleCoefficient a b =
   let bc = Complex.conjugate b
   in  A.lift ((a + bc) / 2, (a - bc) * (-imagUnit / 2))


twoFromSpectrum ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int) (Complex a) ->
   Acc (Array (sh:.Int) (Complex a, Complex a)) ->
   Acc (Array (sh:.Int) (a,a))
twoFromSpectrum subTrans =
   A.map (Exp.modify (expr:+expr) $ \(x:+y) -> (x,y)) .
   subTrans . entangleSpectra

twoFromSpectrum2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Fourier.Transform (sh:.Int:.Int) (Complex a) ->
   Acc (Array (sh:.Int:.Int) (Complex a, Complex a)) ->
   Acc (Array (sh:.Int:.Int) (a,a))
twoFromSpectrum2d subTrans =
   A.map (Exp.modify (expr:+expr) $ \(x:+y) -> (x,y)) .
   subTrans . entangleSpectra2d

entangleSpectra ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int) (Complex a, Complex a)) ->
   Acc (Array (sh:.Int) (Complex a))
entangleSpectra = entangleSpectraGen

entangleSpectra2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array (sh:.Int:.Int) (Complex a, Complex a)) ->
   Acc (Array (sh:.Int:.Int) (Complex a))
entangleSpectra2d = entangleSpectraGen

entangleSpectraGen ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   Acc (Array sh (Complex a, Complex a)) ->
   Acc (Array sh (Complex a))
entangleSpectraGen = A.map (A.fst . A.uncurry entangleCoefficient)


{-
2 *c = a + bc     a  = c + i*d
2i*d = a - bc     bc = c - i*d
-}
entangleCoefficient ::
   (A.RealFloat a, A.Elt (Complex a)) =>
   Exp (Complex a) -> Exp (Complex a) -> Exp (Complex a, Complex a)
entangleCoefficient c d =
   let di = d * imagUnit
   in  A.lift (c + di, Complex.conjugate (c - di))


imagUnit :: (A.Num a, A.Elt (Complex a)) => Exp (Complex a)
imagUnit = A.lift $ zero :+ one

zero, one :: (A.Num a) => Exp a
zero = 0
one = 1
