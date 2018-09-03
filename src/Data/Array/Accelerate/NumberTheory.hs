-- duplicate of synthesizer-core:NumberTheory
module Data.Array.Accelerate.NumberTheory where

import qualified Data.List.HT as ListHT
import Data.Maybe.HT (toMaybe, )
import Data.Bits ((.&.), (.|.), shiftR, )
import Data.List (unfoldr, )


{- |
List all factorizations of an odd number
where the first factor is at most the second factor
and the first factors are in descending order.
-}
fermatFactors :: Integer -> [(Integer,Integer)]
fermatFactors n =
   let root = floor $ sqrt (fromInteger n :: Double)
   in  map (\(a,b) -> (b-a,b+a)) $
       mergeAndFilter
          (zip (scanl (+) n [1,3..]) [0 .. div (n-1) 2])
          (zip (scanl (+) (root*root) $ iterate (2+) (2*root+1)) [root..])

mergeAndFilter :: (Ord a) => [(a,b)] -> [(a,c)] -> [(b,c)]
mergeAndFilter ((a0,b):a0s) ((a1,c):a1s) =
   case compare a0 a1 of
      LT -> mergeAndFilter a0s ((a1,c):a1s)
      GT -> mergeAndFilter ((a0,b):a0s) a1s
      EQ -> (b,c) : mergeAndFilter a0s a1s
mergeAndFilter _ _ = []


multiplicativeGenerator :: Integer -> Integer
multiplicativeGenerator p =
   head $ primitiveRootsOfUnity p (p-1)

primitiveRootsOfUnity :: Integer -> Integer -> [Integer]
primitiveRootsOfUnity modu order =
   let greatDivisors = map (div order) $ uniquePrimeFactors order
   in  filter
          (\n ->
             let pow y = modularPower modu y n
             in  coprime n modu
                 &&
                 pow order == 1
                 &&
                 all (\y -> pow y /= 1) greatDivisors) $
       [1 .. modu-1]


coprime :: Integer -> Integer -> Bool
coprime x y  =  gcd x y == 1

modularPower :: Integer -> Integer -> Integer -> Integer
modularPower modu =
   let go 0 _ = 1
       go expo n =
          case divMod expo 2 of
             (expo2, r) ->
                let n2 = mod (n*n) modu
                in  if r==0
                      then go expo2 n2
                      else mod (go expo2 n2 * n) modu
   in  go

uniquePrimeFactors :: Integer -> [Integer]
uniquePrimeFactors n =
   let oddFactors =
          foldr
             (\p go m ->
                let (q,r) = divMod m p
                in  if r==0
                      then p : go (divideByMaximumPower p q)
                      else
                        if q >= p
                          then go m
                          else if m==1 then [] else m : [])
             (error "uniquePrimeFactors: end of infinite list")
             (iterate (2+) 3)
   in  case powerOfTwoFactors n of
          (1,m) -> oddFactors m
          (_,m) -> 2 : oddFactors m

divideByMaximumPower :: Integer -> Integer -> Integer
divideByMaximumPower b n =
   last $
   n : unfoldr (\m -> case divMod m b of (q,r) -> toMaybe (r==0) (q,q)) n

powerOfTwoFactors :: Integer -> (Integer, Integer)
powerOfTwoFactors n =
   let powerOfTwo = n .&. (-n)
   in  (powerOfTwo, div n powerOfTwo)


ceilingPowerOfTwo :: Integer -> Integer
ceilingPowerOfTwo 0 = 1
ceilingPowerOfTwo n =
   (1+) $ fst $ head $
   dropWhile (uncurry (/=)) $
   ListHT.mapAdjacent (,) $
   scanl (\m d -> shiftR m d .|. m) (n-1) $
   iterate (2*) 1

{-
For every reasonable pair of powers of 3 and 5
it computes the least power of 2,
such that their product is above @n@.
-}
ceiling5Smooth :: Integer -> Integer
ceiling5Smooth n =
   minimum $ map (minimum . ceilingSmooths 2 5 n) $
   ceilingSmooths 2 3 n $ ceilingPowerOfTwo n

{- |
@ceilingSmooths a b n m@
replaces successively @a@ factors in @m@ by @b@ factors
while keeping the product above @n@.
-}
ceilingSmooths :: Integer -> Integer -> Integer -> Integer -> [Integer]
ceilingSmooths a b n =
   let divMany k =
          case divMod k a of
             (q,r) -> if r==0 && q>=n then divMany q else k
       go m  =  m : if mod m a == 0 then go $ divMany $ m*b else []
   in  go
