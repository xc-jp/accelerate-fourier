{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
module Data.Array.Accelerate.Fourier.Utility (
   scaleDown,
   ) where

import           Data.Array.Accelerate.Fourier.Private (Transform)

import           Data.Array.Accelerate.Data.Complex    (Complex ((:+)))

import           Data.Array.Accelerate                 ((:.), Exp, Shape, Slice)
import qualified Data.Array.Accelerate                 as A
import qualified Data.Array.Accelerate.Utility.Sliced  as Sliced


scaleDown ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Transform (sh:.Int) (Complex a)
scaleDown zs =
   A.map (cscale (recip $ A.fromIntegral $ Sliced.length zs)) zs

cscale :: (A.Num a, A.Elt (Complex a)) => Exp a -> Exp (Complex a) -> Exp (Complex a)
cscale x z =
   case A.unlift z of
      r :+ i -> A.lift (x*r :+ x*i)
