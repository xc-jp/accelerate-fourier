{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE TypeOperators    #-}
{- |
The implementations in this module require
that you know the transformation data set size on the Haskell side.
This knowledge is baked into the Accelerate code.
The advantage is,
that you can share preprocessing between calls to the Fourier transforms,
like in:

> let transform = dit2 1024
> in  transform x ... transform y
-}
module Data.Array.Accelerate.Fourier.Preprocessed (
   Transform,
   ditSplitRadix,
   dit2,
   dif2,

   Sign.Sign,
   Sign.forward,
   Sign.inverse,

   transform2d,
   transform3d,

   SubTransformPair(SubTransformPair),
   SubTransformTriple(SubTransformTriple),
   ) where

import           Data.Array.Accelerate.Fourier.Private (PairTransform, SubTransformPair (SubTransformPair),
                                                        SubTransformTriple (SubTransformTriple),
                                                        Transform)
import qualified Data.Array.Accelerate.Fourier.Private as Fourier
import           Data.Array.Accelerate.Fourier.Sign    (Sign)
import qualified Data.Array.Accelerate.Fourier.Sign    as Sign

import           Data.Array.Accelerate.Data.Complex    (Complex)
import           Data.Array.Accelerate.LinearAlgebra   (zipExtrudedVectorWith)
import qualified Data.Array.Accelerate.LinearAlgebra   as LinAlg

import qualified Data.Array.Accelerate.Utility.Sliced  as Sliced

import           Data.Array.Accelerate                 ((:.), Exp, Shape, Slice)
import qualified Data.Array.Accelerate                 as A


{- |
Decimation in time for power-of-two using the split-radix algorithm.
Should be faster than 'dit2'.
-}
ditSplitRadix ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Sign a ->
   Int ->
   Transform (sh:.Int) (Complex a)
ditSplitRadix mode len =
   if len<2
     then id
     else
        Fourier.finishSplitRadix . fst .
        ditSplitRadixGo (A.constant mode) (div len 2) .
        Fourier.initSplitRadix

{- |
Compute the Fourier transforms
of a collection of 2N length signals
and a collection of N length signals
and share some computations between them.
The global extent of @sh@ of all arrays must be equal.
First array must have extent @sh:.count2:.2*len@
and second array must have extent @sh:.count1:.len@.
If this is a restriction for you,
you may use 'Fourier.finishSplitRadixFlat' and 'Fourier.initSplitRadixFlat'
which merge the global shape with our auxiliary dimension
and then work with @sh = Z@.
-}
ditSplitRadixGo ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) ->
   Int ->
   PairTransform (sh:.Int:.Int) (Complex a)
ditSplitRadixGo mode len =
   if len<=1
     then Fourier.ditSplitRadixBase
     else
        let len2 = div len 2
            twiddles = Fourier.twiddleFactorsSRPair mode (A.constant len2)
            imag = Fourier.imagSplitRadix mode
        in  Fourier.ditSplitRadixStep imag twiddles .
            ditSplitRadixGo mode len2 .
            Fourier.ditSplitRadixReorder


{- |
Decimation in time for power-of-two sizes.
-}
dit2 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Sign a ->
   Int ->
   Transform (sh:.Int) (Complex a)
dit2 mode len =
   if len<=1
     then id
     else
        let len2 = div len 2
        in  Fourier.transformRadix2InterleavedTime
               (Fourier.twiddleFactors2 (A.constant mode) (A.constant len2))
               (dit2 mode len2)


{- |
Decimation in frequency for power-of-two sizes.
-}
dif2 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Sign a ->
   Int ->
   Transform (sh:.Int) (Complex a)
dif2 mode len =
   if len<=1
     then id
     else
        let len2 = div len 2
            twiddles = Fourier.twiddleFactors2 (A.constant mode) (A.constant len2)
        in  \arr ->
              let part0 = Sliced.take (A.constant len2) arr
                  part1 = Sliced.drop (A.constant len2) arr
                  evens = A.zipWith (+) part0 part1
                  odds =
                     zipExtrudedVectorWith (*) twiddles $
                     A.zipWith (-) part0 part1
              in  Fourier.merge $ dif2 mode len2 $ Fourier.stack evens odds


{- |
Transforms in 'SubTransformPair'
are ordered from least-significant to most-significant dimension.
-}
transform2d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   SubTransformPair (Complex a) ->
   Transform (sh:.Int:.Int) (Complex a)
transform2d (SubTransformPair transform0 transform1) =
   LinAlg.transpose . transform1 .
   LinAlg.transpose . transform0

{- |
Transforms in 'SubTransformTriple'
are ordered from least-significant to most-significant dimension.
-}
transform3d ::
   (Shape sh, Slice sh, A.RealFloat a, A.Elt (Complex a)) =>
   SubTransformTriple (Complex a) ->
   Transform (sh:.Int:.Int:.Int) (Complex a)
transform3d (SubTransformTriple transform0 transform1 transform2) =
   Fourier.cycleDim3 . transform2 .
   Fourier.cycleDim3 . transform1 .
   Fourier.cycleDim3 . transform0
