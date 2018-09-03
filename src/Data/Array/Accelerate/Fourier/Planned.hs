{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{- |
Like "Data.Array.Accelerate.Fourier.Preprocessed"
this module allows to factor out some preprocessing.
Additionally it gives you concrete objects (plans and caches)
for sharing preprocessed data between transforms.
You cannot only share the preprocessing between transforms of the same size,
but across all array sizes.
This implementation also has the largest collection of algorithms
and thus should be generally fastest among all implementations in this package.
-}
module Data.Array.Accelerate.Fourier.Planned (
   -- * Transforms
   Transform,

   transform,

   transformDecompose,
   transformChirp2,
   transformChirp235,

   convolveCyclic,

   -- * Planning
   Plan,
   plan,
   transformWithPlanner,

   PlanMap,
   planWithMapUpdate,
   planDecomposeWithMapUpdate,
   planChirpWithMapUpdate,
   smallPlanMap,

   -- * Caching
   Cache,
   cache,
   cacheDuplex,
   transformWithCache,

   CacheMap,
   Direction(..),
   cacheFromPlanWithMapUpdate, directionMode,
   cacheFromPlanWithMapUpdate2, directionModes,

   -- * Miscellaneous
   Sign.Sign,
   Sign.forward,
   Sign.inverse,
   ) where

import qualified Data.Array.Accelerate.Convolution.Adhoc as Convolution
import qualified Data.Array.Accelerate.Permutation as Permutation
import qualified Data.Array.Accelerate.NumberTheory as NumberTheory
import qualified Data.Array.Accelerate.Fourier.Private as Fourier
import qualified Data.Array.Accelerate.Fourier.Sign as Sign
import Data.Array.Accelerate.Fourier.Private
          (SubTransform(SubTransform), SubTransformPair(SubTransformPair),
           SubPairTransform(SubPairTransform),
           PairTransform, Transform, )
import Data.Array.Accelerate.Fourier.Utility (scaleDown, )
import Data.Array.Accelerate.Fourier.Sign (Sign(Sign))

import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg
import Data.Array.Accelerate.LinearAlgebra
          (zipExtrudedVectorWith, zipExtrudedMatrixWith, )

import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import Data.Array.Accelerate.Utility.Lift.Exp (expr)

import qualified Data.Array.Accelerate.Utility.Sliced as Sliced
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Data.Complex (Complex, conjugate, )
import Data.Array.Accelerate
          (Exp, Acc, Array, DIM1, DIM2, (:.)((:.)), Elt, Slice, Shape, )

import qualified Control.Monad.Trans.State as State
import Control.Monad (liftM2, )
import Control.Applicative ((<$>), )
import Data.Traversable (for, )

import qualified Data.Map as Map
import Data.Tuple.HT (mapPair, )


{- |
Fourier transform of arbitrary size.
Sign can be

* @forward@: from time domain to frequency spectrum

* @inverse@: from frequency spectrum to time domain

You may share @transform sign n@ between several calls
in order to run some preprocessing only once.
You must make sure that the @length@
is equal to the extent of the inner dimension of every transformed array.
-}
transform ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   Sign a -> Int -> Transform (sh:.Int) (Complex a)
transform sign len = transformWithCache $ cache sign len


{- |
Transform using only Cooley-Tukey, Good-Thomas, Rader, Split-Radix,
but no Bluestein.
This is more for testing and benchmarking than for real use.
-}
transformDecompose ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   Sign a -> Int ->
   Transform (sh :. Int) (Complex a)
transformDecompose =
   transformWithPlanner planDecomposeWithMapUpdate

transformWithPlanner ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   (Integer -> State.State PlanMap Plan) ->
   Sign a -> Int ->
   Transform (sh :. Int) (Complex a)
transformWithPlanner planner sign len =
   transformWithCache $
   cacheFromPlan
      (flip State.evalState smallPlanMap $ planner $ fromIntegral len) $
   directionMode sign len


{- |
The size and type of the signal must match the parameters,
that the cache was generated for.
-}
transformWithCache ::
   (Slice sh, Shape sh, A.RealFloat a) =>
   Cache (Complex a) -> Transform (sh:.Int) (Complex a)
transformWithCache ch =
   case ch of
      CacheIdentity -> id
      CacheSmall size ->
         case size of
            LevelCache2 zs -> Fourier.transform2 zs
            LevelCache3 zs -> Fourier.transform3 zs
            LevelCache4 zs -> Fourier.transform4 zs
            LevelCache5 zs -> Fourier.transform5 zs
      CacheRadix2 level subCache ->
         transformRadix2InterleavedTime level $
         subTransformWithCache subCache
      CacheSplitRadix chain ->
         Fourier.finishSplitRadix . fst .
         transformSplitRadixInterleavedTimeChain chain .
         Fourier.initSplitRadix
      CachePrime level subCaches ->
         transformPrime level $
         fmap subTransformPairWithCache subCaches
      CacheCoprime level subCaches ->
         transformCoprime level $
         subTransformPairWithCache subCaches
      CacheComposite level subCaches ->
         transformComposite level $
         subTransformPairWithCache subCaches
      CacheChirp level subCaches ->
         transformChirp level $
         subTransformPairWithCache subCaches

subTransformWithCache ::
   (A.RealFloat a) => Cache (Complex a) -> SubTransform (Complex a)
subTransformWithCache ch = SubTransform (transformWithCache ch)

subTransformPairWithCache ::
   (A.RealFloat a) =>
   (Cache (Complex a), Cache (Complex a)) -> SubTransformPair (Complex a)
subTransformPairWithCache (ch0,ch1) =
   SubTransformPair (transformWithCache ch0) (transformWithCache ch1)


{- |
Memorize factorizations of the data size and permutation vectors.
-}
data Plan = Plan Integer PlanStructure
   deriving (Show)

data PlanStructure =
     PlanIdentity
   | PlanSmall LevelSmall
   | PlanRadix2 Plan
   | PlanSplitRadix Plan
   | PlanPrime (Maybe Plan)
   | PlanCoprime (Plan, Plan)
   | PlanComposite (Plan, Plan)
   | PlanChirp Plan
   deriving (Show)


{- |
Plan transform algorithms for a certain array size.
-}
plan :: Integer -> Plan
plan n =
   State.evalState (planWithMapUpdate n) smallPlanMap

{- |
Too many nested Rader transformations slow down the transform,
up to quadratic time in the worst case.
As a heuristic we allow at most nesting depth two,
and switch to Bluestein transformation otherwise.
We could compute more precise operation counts
and base our decision on these,
but we found that the actual execution time
differs considerably from the operation counts.
-}
planWithMapUpdate :: Integer -> State.State PlanMap Plan
planWithMapUpdate n = do
   p <- planDecomposeWithMapUpdate n
   if planCountPrimes p < 3
     then return p
     else planChirpWithMapUpdate NumberTheory.ceiling5Smooth n

planCountPrimes :: Plan -> Int
planCountPrimes (Plan _ struct) =
   case struct of
      PlanIdentity -> 0
      PlanSmall _ -> 0
      PlanRadix2 p -> planCountPrimes p
      PlanSplitRadix p -> planCountPrimes p
      PlanPrime mp -> 1 + maybe 0 planCountPrimes mp
      PlanCoprime (m, n) -> max (planCountPrimes m) (planCountPrimes n)
      PlanComposite (m, n) -> max (planCountPrimes m) (planCountPrimes n)
      PlanChirp p -> planCountPrimes p

type PlanMap = Map.Map Integer Plan

{- |
Map of primitive transforms.
You should use this as the initial map
when evaluating a planning sequence using 'State.evalState'.
-}
smallPlanMap :: PlanMap
smallPlanMap =
   Map.fromAscList $ zipWith (\n struct -> (n, Plan n struct)) [0..] $
   PlanIdentity :
   PlanIdentity :
   PlanSmall Level2 :
   PlanSmall Level3 :
   PlanSmall Level4 :
   PlanSmall Level5 :
   []

{- |
Detect and re-use common sub-plans.
-}
planDecomposeWithMap :: Integer -> State.State PlanMap Plan
planDecomposeWithMap n =
   fmap (Plan n) $
   case divMod n 2 of
      (n2,0) ->
         let psr = PlanSplitRadix <$> planDecomposeWithMapUpdate n2
             _pc = PlanComposite <$> planDecomposeWithMapUpdate2 (2,n2)
             _pr2 = PlanRadix2 <$> planDecomposeWithMapUpdate n2
         in  psr
      _ ->
         let facs = NumberTheory.fermatFactors n
         in  -- find unitary divisors
             case filter (\(a,b) -> a>1 && gcd a b == 1) facs of
                q2 : _ -> PlanCoprime <$> planDecomposeWithMapUpdate2 q2
                _ ->
                   let (q2 : _) = facs
                   in  if fst q2 == 1
                         then
                            PlanPrime <$>
                            if False
                              then return Nothing
                              else Just <$> planDecomposeWithMapUpdate (n-1)
                         else PlanComposite <$> planDecomposeWithMapUpdate2 q2

planDecomposeWithMapUpdate :: Integer -> State.State PlanMap Plan
planDecomposeWithMapUpdate n = do
   item <- State.gets (Map.lookup n)
   case item of
      Just p -> return p
      Nothing -> do
         m <- planDecomposeWithMap n
         State.modify (Map.insert n m)
         return m

planDecomposeWithMapUpdate2 ::
   (Integer, Integer) -> State.State PlanMap (Plan, Plan)
planDecomposeWithMapUpdate2 =
   uncurry (liftM2 (,)) .
   mapPair (planDecomposeWithMapUpdate,planDecomposeWithMapUpdate)


{- |
Cache arrays of twiddle factors,
i.e. powers of the primitive root of unity.
-}
data Cache a =
     CacheIdentity
   | CacheSmall (LevelCacheSmall (Exp a))
   | CacheRadix2 (LevelCacheRadix2 a) (Cache a)
   | CacheSplitRadix (CacheSplitRadixChain a)
   | CachePrime (LevelCachePrime a) (Maybe (Cache a, Cache a))
   | CacheCoprime LevelCacheCoprime (Cache a, Cache a)
   | CacheComposite (LevelCacheComposite a) (Cache a, Cache a)
   | CacheChirp (LevelCacheChirp a) (Cache a, Cache a)
   deriving (Show)

data CacheSplitRadixChain a =
     CacheSplitRadixCons (LevelCacheSplitRadix a) (CacheSplitRadixChain a)
   | CacheSplitRadixEnd (Cache a) (Cache a)
   deriving (Show)

{- |
The expression @cache sign len@
precomputes all data that is needed for Fourier transforms
for signals of length @len@.
You can use this cache in 'transformWithCache'.
-}
cache ::
   (A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   Sign a -> Int -> Cache (Complex a)
cache sign len =
   cacheFromPlan
      (plan $ fromIntegral len)
      (directionMode sign len)


{- |
It is @(cache inverse x, cache forward x) = cacheDuplex x@
but 'cacheDuplex' shares common data of both caches.
-}
cacheDuplex ::
   (a ~ Complex b, A.RealFloat b, A.FromIntegral Int b, Num b) =>
   Int -> (Cache a, Cache a)
cacheDuplex len =
   let p = plan $ fromIntegral len
   in  flip State.evalState Map.empty $
          cacheFromPlanWithMapUpdate2 (p,p) (directionModes len)


data Direction = Forward | Inverse
   deriving (Show, Eq, Ord)


type CacheMap a = Map.Map (Integer,Direction) (Cache a)

cacheFromPlan ::
   (a ~ Complex b, A.RealFloat b, A.FromIntegral Int b, Num b) =>
   Plan -> (Direction, Sign b) -> Cache a
cacheFromPlan p z =
   State.evalState (cacheFromPlanWithMapUpdate p z) Map.empty


{- |
Detect and re-use common sub-caches.
-}
cacheFromPlanWithMap ::
   (a ~ Complex b, A.RealFloat b, Num b, A.FromIntegral Int b) =>
   Plan -> (Direction, Sign b) ->
   State.State (CacheMap a) (Cache a)
cacheFromPlanWithMap (Plan len struct) dsign@(_d,sign) =
   case struct of
      PlanIdentity -> return $ CacheIdentity
      PlanSmall size -> return $ CacheSmall $
         case size of
            Level2 -> LevelCache2 $ Fourier.cache2 $ A.constant sign
            Level3 -> LevelCache3 $ Fourier.cache3 $ A.constant sign
            Level4 -> LevelCache4 $ Fourier.cache4 $ A.constant sign
            Level5 -> LevelCache5 $ Fourier.cache5 $ A.constant sign
      PlanRadix2 subPlan@(Plan len2 _) ->
         CacheRadix2 (levelCacheRadix2 len2 sign) <$>
         cacheFromPlanWithMapUpdate subPlan dsign
      PlanSplitRadix subPlan@(Plan len2 subStruct) -> do
         subCache <- cacheFromPlanWithMapUpdate subPlan dsign
         case subCache of
            CacheSplitRadix chain ->
               return $
                  CacheSplitRadix $
                  CacheSplitRadixCons (levelCacheSplitRadix len2 sign) chain
            _ ->
               case subStruct of
                  PlanSplitRadix subsubPlan -> do
                     subsubCache <- cacheFromPlanWithMapUpdate subsubPlan dsign
                     return $ CacheSplitRadix $
                        CacheSplitRadixCons (levelCacheSplitRadix len2 sign) $
                        CacheSplitRadixEnd subCache subsubCache
                  _ ->
                     return $ CacheRadix2 (levelCacheRadix2 len2 sign) subCache
      PlanPrime maybeSubPlan ->
         (\maybeSubCaches ->
            CachePrime
               (levelCachePrime len
                  (fmap (subTransformWithCache . fst) maybeSubCaches) sign)
               maybeSubCaches)
         <$>
         for maybeSubPlan
            (\subPlan ->
               cacheFromPlanWithMapUpdate2 (subPlan,subPlan)
                  (directionModes $ fromInteger $ len-1))
      PlanCoprime subPlans@(Plan n _, Plan m _) ->
         CacheCoprime (levelCacheCoprime (n,m)) <$>
         cacheFromPlanWithMapUpdate2 subPlans (dsign, dsign)
      PlanComposite subPlans@(Plan n _, Plan m _) ->
         CacheComposite (levelCacheComposite (n,m) sign)
         <$>
         cacheFromPlanWithMapUpdate2 subPlans (dsign, dsign)
      PlanChirp subPlan@(Plan padlen _) ->
         (\subCaches ->
            CacheChirp
               (levelCacheChirp len padlen
                  (subTransformWithCache (fst subCaches)) sign)
               subCaches)
         <$>
         cacheFromPlanWithMapUpdate2 (subPlan,subPlan)
            (directionModes $ fromInteger padlen)

cacheFromPlanWithMapUpdate ::
   (a ~ Complex b, A.RealFloat b, A.FromIntegral Int b, Num b) =>
   Plan -> (Direction, Sign b) ->
   State.State (CacheMap a) (Cache a)
cacheFromPlanWithMapUpdate p@(Plan len _) z = do
   let key = (len, fst z)
   item <- State.gets (Map.lookup key)
   case item of
      Just c -> return c
      Nothing -> do
         m <- cacheFromPlanWithMap p z
         State.modify (Map.insert key m)
         return m

cacheFromPlanWithMapUpdate2 ::
   (a ~ Complex b, A.RealFloat b, A.FromIntegral Int b, Num b) =>
   (Plan, Plan) -> ((Direction, Sign b), (Direction, Sign b)) ->
   State.State (CacheMap a) (Cache a, Cache a)
cacheFromPlanWithMapUpdate2 (p0,p1) (dm0,dm1) =
   liftM2 (,)
      (cacheFromPlanWithMapUpdate p0 dm0)
      (cacheFromPlanWithMapUpdate p1 dm1)


directionMode ::
   (Num a, Ord a) =>
   Sign a -> Int -> (Direction, Sign a)
directionMode (Sign sign) len =
   (if sign>0 then fst else snd) $ directionModes len

directionModes ::
   (Num a) =>
   Int -> ((Direction, Sign a), (Direction, Sign a))
directionModes _len =
   ((Inverse, Sign.inverse), (Forward, Sign.forward))



data LevelSmall = Level2 | Level3 | Level4 | Level5
   deriving (Show, Eq, Ord, Enum)

data LevelCacheSmall a =
     LevelCache2 a
   | LevelCache3 (a,a)
   | LevelCache4 (a,a,a)
   | LevelCache5 (a,a,a,a)
   deriving (Show)



data LevelCacheRadix2 a = LevelCacheRadix2 (Acc (Array DIM1 a))
   deriving (Show)

levelCacheRadix2 ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Integer -> Sign a -> LevelCacheRadix2 (Complex a)
levelCacheRadix2 n2 sign =
   LevelCacheRadix2 $
   Fourier.twiddleFactors2 (A.constant sign) (expInteger n2)

transformRadix2InterleavedTime ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCacheRadix2 a ->
   SubTransform a ->
   Transform (sh:.Int) a
transformRadix2InterleavedTime
      (LevelCacheRadix2 twiddles) (SubTransform subTrans) =
   Fourier.transformRadix2InterleavedTime twiddles subTrans


data LevelCacheSplitRadix a =
   LevelCacheSplitRadix a (Acc (Array DIM1 a), Acc (Array DIM1 a))
   deriving (Show)

levelCacheSplitRadix ::
   (A.RealFloat a, Num a, A.FromIntegral Int a) =>
   Integer -> Sign a -> LevelCacheSplitRadix (Complex a)
levelCacheSplitRadix n2 sign =
   LevelCacheSplitRadix (Fourier.imagSplitRadixPlain sign) $
   Fourier.twiddleFactorsSRPair (A.constant sign) (expInteger (div n2 2))

transformSplitRadixInterleavedTime ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCacheSplitRadix a ->
   SubPairTransform a ->
   PairTransform (sh:.Int:.Int) a
transformSplitRadixInterleavedTime
      (LevelCacheSplitRadix imag twiddles) (SubPairTransform subTrans) =
   Fourier.ditSplitRadixStep (A.constant imag) twiddles .
   subTrans .
   Fourier.ditSplitRadixReorder

transformSplitRadixInterleavedTimeChain ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   CacheSplitRadixChain a ->
   PairTransform (sh:.Int:.Int) a
transformSplitRadixInterleavedTimeChain chain =
   case chain of
      CacheSplitRadixCons level remChain ->
         transformSplitRadixInterleavedTime level $
         SubPairTransform (transformSplitRadixInterleavedTimeChain remChain)
      CacheSplitRadixEnd subCache2 subCache1 ->
         mapPair (transformWithCache subCache2, transformWithCache subCache1)


newtype LevelCacheComposite a =
   LevelCacheComposite (Acc (Array DIM2 a))
   deriving (Show)

levelCacheComposite ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   (Integer, Integer) -> Sign a -> LevelCacheComposite (Complex a)
levelCacheComposite (n,m) sign =
   LevelCacheComposite $
   Fourier.twiddleFactors (A.constant sign) (expInteger n) (expInteger m)


{- |
For @transformComposite z (n,m) sig@,
it must hold @n*m == length sig@ and @z ^ length sig == 1@.

Cooley-Tukey-algorithm
-}
transformComposite ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCacheComposite a ->
   SubTransformPair a ->
   Transform (sh:.Int) a
transformComposite
      (LevelCacheComposite twiddles)
      (SubTransformPair subTransN subTransM) =

   Fourier.merge .
   subTransN .
   LinAlg.transpose .
   zipExtrudedMatrixWith (*) twiddles .
   subTransM .
   Sliced.sliceHorizontal (A.shape twiddles)


newtype LevelCacheCoprime = LevelCacheCoprime (Integer, Integer)
   deriving (Show)

{-
Fourier exponent matrix of a signal of size 6.

0 0 0 0 0 0     0               0   0   0       0     0
0 1 2 3 4 5               0       2   0   4       3     0
0 2 4 0 2 4  =          0    *  0   2   4    *      0     0
0 3 0 3 0 3           0           0   0   0     0     3
0 4 2 0 4 2         0           0   4   2         0     0
0 5 4 3 2 1       0               4   0   2         0     3
-}
levelCacheCoprime :: (Integer, Integer) -> LevelCacheCoprime
levelCacheCoprime = LevelCacheCoprime


{- |
For @transformCoprime z (n,m) sig@,
the parameters @n@ and @m@ must be relatively prime
and @n*m == length sig@ and @z ^ length sig == 1@.

Good-Thomas algorithm
-}
transformCoprime ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCacheCoprime ->
   SubTransformPair a ->
   Transform (sh:.Int) a
transformCoprime
      (LevelCacheCoprime (n,m)) (SubTransformPair subTransN subTransM) =
   permuteSkewGridInv .
   subTransM .
   LinAlg.transpose .
   subTransN .
   permuteSkewGrid (expInteger m) (expInteger n)

permuteSkewGrid ::
   (Slice sh, Shape sh, Elt a) =>
   Exp Int -> Exp Int -> LinAlg.Vector sh a -> LinAlg.Matrix sh a
permuteSkewGrid m n arr =
   let (sh:.nm) = Exp.unlift (expr:.expr) $ A.shape arr
   in  A.backpermute
          (A.lift (sh :. m :. n))
          (Exp.modify (expr:.expr:.expr) $
           \(ix:.k:.j) -> ix :. mod (n*k + m*j) nm)
          arr

permuteSkewGridInv ::
   (Slice sh, Shape sh, Elt a) =>
   LinAlg.Matrix sh a -> LinAlg.Vector sh a
permuteSkewGridInv arr =
   let (sh:.m:.n) = Exp.unlift (expr:.expr:.expr) $ A.shape arr
   in  A.backpermute
          (A.lift (sh :. n*m))
          (Exp.modify (expr:.expr) $
           \(ix:.k) -> ix :. mod k m :. mod k n)
          arr



{-
Fourier exponent matrix of a signal of size 7.

0 0 0 0 0 0 0
0 1 2 3 4 5 6
0 2 4 6 1 3 5
0 3 6 2 5 1 4
0 4 1 5 2 6 3
0 5 3 1 6 4 2
0 6 5 4 3 2 1

multiplicative generator in Z7: 3
permutation of rows and columns by powers of 3: 1 3 2 6 4 5

0 0 0 0 0 0 0
0 1 3 2 6 4 5
0 3 2 6 4 5 1
0 2 6 4 5 1 3
0 6 4 5 1 3 2
0 4 5 1 3 2 6
0 5 1 3 2 6 4

Inverse permutation: 1 3 2 5 6 4
The inverse permutations seems not to be generated by a multiplication.
-}
data LevelCachePrime a =
   LevelCachePrime (Permutation.T, Permutation.T) (Acc (Array DIM1 a))
      deriving (Show)

levelCachePrime ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Integer ->
   Maybe (SubTransform (Complex a)) ->
   Sign a -> LevelCachePrime (Complex a)
levelCachePrime n maybeSubTrans sign =
   let len = fromInteger n
       perm = A.use $ Permutation.multiplicative len
       kernel =
          A.map (Sign.cisRat (A.constant sign) (A.constant len) . (1+)) perm
   in  LevelCachePrime
          (Permutation.reverse perm, Permutation.inverse perm)
          (maybe id
             (\(SubTransform subTrans) -> scaleDown . subTrans)
             maybeSubTrans kernel)

{- |
Rader's algorithm for prime length signals.
-}
transformPrime ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCachePrime a ->
   Maybe (SubTransformPair a) ->
   Transform (sh:.Int) a
transformPrime (LevelCachePrime (rev, inv) zs) maybeSubTranss =
   let conv =
          case maybeSubTranss of
             Nothing ->
                \xs ->
                   Convolution.complex
                      (Convolution.cyclic Convolution.karatsuba)
                      (LinAlg.extrudeVector (A.indexTail $ A.shape xs) zs)
                      xs
             Just subTranss ->
                convolveSingleSpectrumCyclicCache subTranss zs
   in  \arr ->
          let x0 = Sliced.head arr
              res = Sliced.tail arr
          in  LinAlg.zipScalarVectorWith (+) x0 $
              Sliced.cons (A.fold1 (+) res) $
              Permutation.apply inv $
              conv $
              Permutation.apply rev res



{- |
Fourier transform for arbitrary lengths
based on the Bluestein transform or chirp z-transform
on an array with power-of-two size.
It may be faster than 'transform' for certain prime factors.
Find bad factors e.g. in <http://oeis.org/A061092> and <http://oeis.org/A059411>
and nicer factors in <http://oeis.org/A061303>.
-}
transformChirp2 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   Sign a -> Int ->
   Transform (sh :. Int) (Complex a)
transformChirp2 = transformChirpComplete NumberTheory.ceilingPowerOfTwo

{- |
Fourier transform for arbitrary lengths
based on the Bluestein transform
on an array with 5-smooth size.
(5-smooth = all prime factors are at most 5)
-}
transformChirp235 ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   Sign a -> Int ->
   Transform (sh :. Int) (Complex a)
transformChirp235 = transformChirpComplete NumberTheory.ceiling5Smooth


transformChirpComplete ::
   (Slice sh, Shape sh, A.RealFloat a, A.FromIntegral Int a, Num a, Ord a) =>
   (Integer -> Integer) ->
   Sign a -> Int ->
   Transform (sh :. Int) (Complex a)
transformChirpComplete padLength =
   transformWithPlanner (planChirpWithMapUpdate padLength)

planChirpWithMapUpdate ::
   (Integer -> Integer) -> Integer -> State.State PlanMap Plan
planChirpWithMapUpdate padLength len =
   Plan len <$>
   if len<2
     then return PlanIdentity
     else PlanChirp <$> planDecomposeWithMapUpdate (padLength (2*len-1))


data LevelCacheChirp a =
   LevelCacheChirp (Acc (Array DIM1 a)) (Acc (Array DIM1 a))
      deriving (Show)

levelCacheChirp ::
   (A.RealFloat a, A.FromIntegral Int a) =>
   Integer -> Integer ->
   SubTransform (Complex a) ->
   Sign a -> LevelCacheChirp (Complex a)
levelCacheChirp len padlen (SubTransform subTrans) sign =
   let chirp =
          Fourier.chirp (A.constant sign) (expInteger padlen) (expInteger len)
   in  LevelCacheChirp
          (A.take (expInteger len) chirp)
          (scaleDown $ subTrans $ A.map conjugate chirp)

{- |
Bluestein's algorithm for signals of arbitrary length
and possibly slightly generalised basis vectors.
-}
transformChirp ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   LevelCacheChirp a ->
   SubTransformPair a ->
   Transform (sh:.Int) a
transformChirp (LevelCacheChirp chirp chirpSpec) subTranss =
   let conv = convolveSingleSpectrumCyclicCache subTranss chirpSpec
       twistChirp = zipExtrudedVectorWith (*) chirp
   in  \arr ->
          twistChirp $
          Sliced.take (Sliced.length arr) $
          conv $
          Sliced.pad 0 (A.length chirpSpec) $
          twistChirp arr


{- |
Signals must have equal size and must not be empty.
-}
convolveCyclic ::
   (Shape sh, Slice sh, a ~ Complex b,
    A.RealFloat b, A.FromIntegral Int b, Num b) =>
   Int ->
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int) a)
convolveCyclic leni =
   let len = fromIntegral leni
       (z,zInv) = directionModes leni
   in  convolveCyclicCache $
       subTransformPairWithCache
          (cacheFromPlan (plan len) z,
           cacheFromPlan (plan len) zInv)

convolveCyclicCache ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b, A.FromIntegral Int b) =>
   SubTransformPair a ->
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int) a) ->
   Acc (Array (sh:.Int) a)
convolveCyclicCache transs@(SubTransformPair trans _) x =
   convolveSpectrumCyclicCache transs $ scaleDown $ trans x

convolveSingleSpectrumCyclicCache ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   SubTransformPair a ->
   Acc (Array DIM1 a) -> Transform (sh:.Int) a
convolveSingleSpectrumCyclicCache caches x y =
   convolveSpectrumCyclicCache caches
      (LinAlg.extrudeVector (A.indexTail $ A.shape y) x) y

{- |
This function does not apply scaling.
That is you have to scale the spectrum by @recip (length x)@
if you want a plain convolution.
-}
convolveSpectrumCyclicCache ::
   (Shape sh, Slice sh, a ~ Complex b, A.RealFloat b) =>
   SubTransformPair a ->
   Acc (Array (sh:.Int) a) -> Transform (sh:.Int) a
convolveSpectrumCyclicCache (SubTransformPair trans transInv) x y =
   transInv $ A.zipWith (*) x (trans y)


expInteger :: (A.Num a) => Integer -> Exp a
expInteger = A.fromInteger
