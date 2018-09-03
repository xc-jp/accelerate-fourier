{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.Fourier.Utility (
   scaleDown,
   ) where

import Data.Array.Accelerate.Fourier.Private (Transform, )

import Data.Array.Accelerate.Data.Complex (Complex((:+)), )

import qualified Data.Array.Accelerate.Utility.Sliced as Sliced
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Exp, Slice, Shape, (:.), )


scaleDown ::
   (Shape sh, Slice sh, A.RealFloat a, A.FromIntegral Int a) =>
   Transform (sh:.Int) (Complex a)
scaleDown zs =
   A.map (cscale (recip $ A.fromIntegral $ Sliced.length zs)) zs

cscale :: (A.Num a) => Exp a -> Exp (Complex a) -> Exp (Complex a)
cscale x z =
   case A.unlift z of
      r :+ i -> A.lift (x*r :+ x*i)
