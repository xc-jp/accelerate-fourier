{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
module Data.Array.Accelerate.Convolution.Preprocessed (
   Transform2,
   karatsuba,
   ) where

import           Data.Array.Accelerate.Convolution.Private (Transform2,
                                                            indexPad)

import           Data.Array.Accelerate.Utility.Lift.Exp    (expr)
import qualified Data.Array.Accelerate.Utility.Lift.Exp    as Exp
import qualified Data.Array.Accelerate.Utility.Sliced      as Sliced

import           Data.Array.Accelerate                     ((:.) ((:.)),
                                                            All (All),
                                                            Any (Any), Shape,
                                                            Slice)
import qualified Data.Array.Accelerate                     as A


{- |
Both arrays must have the same size.

There is not much to preprocess,
thus you should prefer 'Data.Array.Accelerate.Convolution.Adhoc.karatsuba'.
-}
karatsuba ::
   (Shape sh, Slice sh, A.Num a) =>
   Int -> Transform2 (sh :. Int) a
karatsuba len x y =
   if len <= 1
     then A.zipWith (*) x y
     else
        let len2 = - div (-len) 2
            elen2 = A.constant len2
            xl = Sliced.take elen2 x
            yl = Sliced.take elen2 y
            xr = Sliced.pad 0 elen2 $ Sliced.drop elen2 x
            yr = Sliced.pad 0 elen2 $ Sliced.drop elen2 y
            zmerged =
               karatsuba len2
                  (Sliced.stack3 xl (A.zipWith (+) xl xr) xr)
                  (Sliced.stack3 yl (A.zipWith (+) yl yr) yr)
            zl = A.slice zmerged $ A.lift $ Any :. (0::Int) :. All
            zm = A.slice zmerged $ A.lift $ Any :. (1::Int) :. All
            zr = A.slice zmerged $ A.lift $ Any :. (2::Int) :. All
            zc = A.zipWith (-) zm $ A.zipWith (+) zl zr
            sh = A.indexTail $ A.shape zc
        in  A.generate (A.lift $ sh :. 2*len-1) $
            Exp.modify (expr:.expr) $
            \(ix:.k) ->
               indexPad (ix:.k)         zl +
               indexPad (ix:.k-elen2)   zc +
               indexPad (ix:.k-elen2*2) zr
