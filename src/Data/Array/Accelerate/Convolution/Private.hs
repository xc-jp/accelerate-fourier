{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.Convolution.Private where

import qualified Data.Array.Accelerate.Utility.Sliced as Sliced

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Exp, Acc, Array, (:.)((:.)), Slice, Shape, (!), (?), )

import Prelude (Int, )


type Transform2 sh a =
   Acc (Array sh a) ->
   Acc (Array sh a) ->
   Acc (Array sh a)


indexPad ::
   (Shape sh, Slice sh, A.Num a) =>
   Exp sh :. Exp Int ->
   Acc (Array (sh:.Int) a) -> Exp a
indexPad (ix:.k) xs =
   0 A.<= k A.&& k A.< Sliced.length xs ? (xs ! A.lift (ix:.k), 0)
