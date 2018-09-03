{- |
Permutations of signals as needed for Fast Fourier transforms.
Most functions are independent of the Signal framework.
We could move them as well to Synthesizer.Basic.
-}
module Data.Array.Accelerate.Permutation where

import qualified Data.Array.Accelerate.NumberTheory as NumberTheory
import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg

import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import Data.Array.Accelerate.Utility.Lift.Exp (expr)

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Exp, Acc, Array, DIM1, Elt, Z(Z), (:.)((:.)),
           Slice, Shape, (!), )


type T = Acc Plain
type Plain = Array DIM1 Int


apply ::
   (Slice sh, Shape sh, Elt a) =>
   T -> LinAlg.Vector sh a -> LinAlg.Vector sh a
apply p xs =
   A.generate
      (Exp.modify2 (expr:.expr) (expr:.expr)
         (\(_z:.n) (sh:._) -> (sh:.n))
         (A.shape p) (A.shape xs)) $
   Exp.modify (expr:.expr) $
   \(ix:.k) -> xs ! A.lift (ix :. p ! A.index1 k)


plainSize :: Plain -> Int
plainSize arr =
   case A.arrayShape arr of
      _ :. n -> n

size :: T -> Exp Int
size = A.length


{- |
Beware of 0-based indices stored in the result vector.
-}
multiplicative :: Int -> Plain
multiplicative ni =
   let n = fromIntegral ni
       gen = NumberTheory.multiplicativeGenerator n
   in  A.fromList (Z :. ni-1) $
       map (fromInteger . subtract 1) $
       iterate (\x -> mod (gen * x) n) 1


{- |
We only need to compute the inverse permutation explicitly,
because not all signal structures support write to arbitrary indices,
thus Generic.Write does not support it.
For strict StorableVector it would be more efficient
to build the vector directly.

It holds:

> inverse . inverse == id
-}
inverse :: T -> T
inverse perm =
   A.permute (+)
      (A.fill (A.shape perm) 0)
      (A.index1 . (perm!))
      (A.generate (A.shape perm) A.unindex1)

reverse :: T -> T
reverse perm =
   A.backpermute
      (A.shape perm)
      (\ix -> A.index1 $ mod (- A.unindex1 ix) (A.unindex1 $ A.shape perm))
      perm
