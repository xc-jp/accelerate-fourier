module Main where

import Criterion.Main (Benchmark, defaultMain, bgroup, bench, whnf, )

import qualified Data.Array.Accelerate.Fourier.Planned as Planned
import qualified Data.Array.Accelerate.Fourier.Preprocessed as Prep
import qualified Data.Array.Accelerate.Fourier.Adhoc as Adhoc

import qualified Data.Array.Accelerate.LLVM.Native as LLVM
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate (Array, DIM2, Z(Z), (:.)((:.)), )
import Data.Complex (Complex, )


powersOfTwo ::
   (Int ->
    A.Acc (Array DIM2 (Complex Float)) ->
    A.Acc (Array DIM2 (Complex Float))) ->
   [Benchmark]
powersOfTwo f =
   take 6 $ flip map (iterate (2*) 1024) $ \len ->
      bench (show len) $ whnf (LLVM.run1 (f len)) $
      A.fromList (Z:.16:.len) $ repeat 0

powersOfTwos :: [Benchmark]
powersOfTwos =
   bgroup "split-radix adhoc"
      (powersOfTwo (const $ Adhoc.ditSplitRadix Adhoc.forward)) :
   bgroup "split-radix preprocessed"
      (powersOfTwo (Prep.ditSplitRadix Prep.forward)) :
   bgroup "dit2 adhoc"
      (powersOfTwo (const $ Adhoc.dit2 Adhoc.forward)) :
   bgroup "dit2 preprocessed"
      (powersOfTwo (Prep.dit2 Prep.forward)) :
   bgroup "dif2 preprocessed"
      (powersOfTwo (Prep.dif2 Prep.forward)) :
   bgroup "dit2 planned"
      (powersOfTwo (Planned.transform Planned.forward)) :
   []


arbitrary ::
   (Int ->
    A.Acc (Array DIM2 (Complex Float)) ->
    A.Acc (Array DIM2 (Complex Float))) ->
   [Benchmark]
arbitrary f =
   take 128 $ flip map (iterate (1+) 1) $ \len ->
      bench (show len) $ whnf (LLVM.run1 (f len)) $
      A.fromList (Z:.4096:.len) $ repeat 0

arbitraryLengths :: [Benchmark]
arbitraryLengths =
   bgroup "auto planned"
      (arbitrary (Planned.transform Planned.forward)) :
   bgroup "decompose planned"
      (arbitrary (Planned.transformDecompose Planned.forward)) :
   bgroup "chirp235 planned"
      (arbitrary (Planned.transformChirp235 Planned.forward)) :
   bgroup "chirp2 planned"
      (arbitrary (Planned.transformChirp2 Planned.forward)) :
   bgroup "auto adhoc"
      (arbitrary (const $ Adhoc.transform Adhoc.forward)) :
   bgroup "chirp235 adhoc"
      (arbitrary (const $ Adhoc.transformChirp235 Adhoc.forward)) :
   bgroup "chirp2 adhoc"
      (arbitrary (const $ Adhoc.transformChirp2 Adhoc.forward)) :
   []


main :: IO ()
main = defaultMain $
   bgroup "2^n" powersOfTwos :
   bgroup "any" arbitraryLengths :
   []
