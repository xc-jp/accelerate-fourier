{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Data.Array.Accelerate.Cyclic (
   Transform,
   reverse,
   reverse2d,
   ) where

import Data.Array.Accelerate.Fourier.Private (Transform)

import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import Data.Array.Accelerate.Utility.Lift.Exp (expr)

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate ((:.)((:.)), )

import Prelude (Int, mod, ($), )


reverse ::
   (A.Shape sh, A.Slice sh, A.Elt a) =>
   Transform (sh :. Int) a
reverse arr =
   let sh = A.shape arr
   in  A.backpermute sh
          (Exp.modify (expr:.expr) $
           \(ix:.k) -> ix :. mod (-k) (A.indexHead sh))
          arr

reverse2d ::
   (A.Shape sh, A.Slice sh, A.Elt a) =>
   Transform (sh :. Int :. Int) a
reverse2d arr =
   let sh = A.shape arr
       (_z:.height:.width) = Exp.unlift (expr:.expr:.expr) sh
   in  A.backpermute sh
          (Exp.modify (expr:.expr:.expr) $
           \(ix:.y:.x) -> ix :. mod (-y) height :. mod (-x) width)
          arr
