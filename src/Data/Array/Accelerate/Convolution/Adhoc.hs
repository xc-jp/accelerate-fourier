{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
module Data.Array.Accelerate.Convolution.Adhoc (
   Transform2,
   karatsuba,
   cyclic,
   complex,
   ) where

import           Data.Array.Accelerate.Convolution.Private (Transform2,
                                                            indexPad)
import           Data.Array.Accelerate.Fourier.Private     (Transform)

import           Data.Array.Accelerate.Utility.Lift.Acc    (acc)
import qualified Data.Array.Accelerate.Utility.Lift.Acc    as Acc
import           Data.Array.Accelerate.Utility.Lift.Exp    (expr)
import qualified Data.Array.Accelerate.Utility.Lift.Exp    as Exp
import qualified Data.Array.Accelerate.Utility.Sliced      as Sliced
import qualified Data.Array.Accelerate.Utility.Sliced1     as Sliced1

import           Data.Array.Accelerate.Data.Complex        (Complex ((:+)))
import qualified Data.Array.Accelerate.Data.Complex        as Complex

import           Data.Array.Accelerate                     ((:.) ((:.)), Acc,
                                                            All (All),
                                                            Any (Any), Array,
                                                            Exp, Shape, Slice,
                                                            Z (Z), (!))
import qualified Data.Array.Accelerate                     as A


{- |
Both arrays must have the same size.
-}
karatsuba ::
   (Shape sh, Slice sh, A.Num a) =>
   Transform2 (sh :. Int) a
karatsuba x y =
   flip A.slice (A.lift $ Any :. (0::Int) :. All)
   .
   A.afst
   .
   A.awhile
      (\arrs -> A.unit $ (Sliced.length $ A.asnd arrs) A.> 1)
      (Acc.modify (acc, acc) $
       \(z, lens) ->
          (karatsubaGo (lens ! A.index1 0) (2*(lens ! A.index1 1)-1) z,
           Sliced.tail lens))
   .
   (Acc.modify ((acc, acc), acc) $
    \((x0,y0), lens) -> (A.zipWith (*) x0 y0, lens))
   .
   A.awhile
      (\arrs -> A.unit $ (Sliced.length $ A.afst $ A.afst arrs) A.> 1)
      (Acc.modify ((acc, acc), acc) $
       \((x0,y0), lens) ->
          let (x1,y1) = karatsubaReorder (x0,y0)
          in  ((x1,y1), Sliced.consExp (Sliced.length x1) lens))
   $
   A.lift
      ((A.replicate (A.lift $ Any :. (1::Int) :. All) x,
        A.replicate (A.lift $ Any :. (1::Int) :. All) y),
       A.fill (A.constant $ Z:.1) (Sliced.length x))

karatsubaReorder ::
   (Shape sh, Slice sh, A.Num a) =>
   (Acc (Array (sh :. Int :. Int) a),
    Acc (Array (sh :. Int :. Int) a)) ->
   (Acc (Array (sh :. Int :. Int) a),
    Acc (Array (sh :. Int :. Int) a))
karatsubaReorder (x,y) =
   let len2 = - div (- Sliced.length x) 2
       xl = Sliced.take len2 x
       yl = Sliced.take len2 y
       xr = Sliced.pad 0 len2 $ Sliced.drop len2 x
       yr = Sliced.pad 0 len2 $ Sliced.drop len2 y
   in  (Sliced1.append3 xl (A.zipWith (+) xl xr) xr,
        Sliced1.append3 yl (A.zipWith (+) yl yr) yr)

karatsubaGo ::
   (Shape sh, Slice sh, A.Num a) =>
   Exp Int ->
   Exp Int ->
   Transform (sh :. Int :. Int) a
karatsubaGo xlen zlen zmerged =
   let (sh:.n:._m) = Exp.unlift (expr:.expr:.expr) $ A.shape zmerged
       n3 = div n 3
       zl = Sliced1.take n3 zmerged
       zm = Sliced1.drop n3 zmerged
       zr = Sliced1.drop (2*n3) zmerged
       zc = A.zipWith (-) zm $ A.zipWith (+) zl zr
   in  A.generate (A.lift $ sh :. n3 :. zlen) $
       Exp.modify (expr:.expr) $
       \(ix:.k) ->
          indexPad (ix:.k)        zl +
          indexPad (ix:.k-xlen)   zc +
          indexPad (ix:.k-xlen*2) zr


{- |
Turn an ordinary convolution into a cyclic convolution of the same length.
-}
cyclic ::
   (Shape sh, Slice sh, A.Num a) =>
   Transform2 (sh :. Int) a ->
   Transform2 (sh :. Int) a
cyclic conv x y =
   let z = conv x y
       len = Sliced.length x
   in  A.zipWith (+) z $ Sliced.pad 0 len $ Sliced.drop len z


{- |
Turn a real-valued convolution into a complex-valued convolution.
Can be removed when we get @instance IsNum (Complex a)@.
-}
complex, _complex ::
   (Shape sh, Slice sh, A.Num a, A.Elt (Complex a)) =>
   Transform2 (sh :. Int) a ->
   Transform2 (sh :. Int) (Complex a)
complex conv x y =
   let xr = A.map Complex.real x; xi = A.map Complex.imag x
       yr = A.map Complex.real y; yi = A.map Complex.imag y
       xm = A.zipWith (+) xr xi
       ym = A.zipWith (+) yr yi
       xryr = conv xr yr
       xiyi = conv xi yi
       xmym = conv xm ym
   in  A.zipWith
          (Exp.modify2 expr expr (:+))
          (A.zipWith (-) xryr xiyi)
          (A.zipWith (-) xmym $ A.zipWith (+) xryr xiyi)

_complex conv x y =
   let xr = A.map Complex.real x; xi = A.map Complex.imag x
       yr = A.map Complex.real y; yi = A.map Complex.imag y
   in  A.zipWith
          (Exp.modify2 expr expr (:+))
          (A.zipWith (-) (conv xr yr) (conv xi yi))
          (A.zipWith (+) (conv xr yi) (conv xi yr))
