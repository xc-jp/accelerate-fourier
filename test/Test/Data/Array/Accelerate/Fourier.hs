{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Test.Data.Array.Accelerate.Fourier where -- (tests) where

import qualified Data.Array.Accelerate.Fourier.Real as FourierReal
import qualified Data.Array.Accelerate.Fourier.Preprocessed as Prep
import qualified Data.Array.Accelerate.Fourier.Adhoc as Adhoc
import qualified Data.Array.Accelerate.Fourier.Planned as Planned
import qualified Data.Array.Accelerate.Convolution.Preprocessed as ConvPrep
import qualified Data.Array.Accelerate.Convolution.Adhoc as Convolution
import qualified Data.Array.Accelerate.Cyclic as Cyclic
import qualified Data.Array.Accelerate.Interpreter as AI
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Fourier.Planned (Transform, )
import Data.Array.Accelerate.Fourier.Utility (scaleDown, )
import Data.Array.Accelerate
          (Acc, Exp, Array, DIM1, DIM2, DIM3, Z(Z), (:.)((:.)), )

import qualified Data.Array.Accelerate.LinearAlgebra as LinAlg
import qualified Data.Array.Accelerate.Utility.Sliced as Sliced
import qualified Data.Array.Accelerate.Utility.Lift.Acc as Acc
import qualified Data.Array.Accelerate.Utility.Lift.Exp as Exp
import Data.Array.Accelerate.Utility.Lift.Exp (expr)

import qualified Test.QuickCheck as QC
import Test.QuickCheck (Arbitrary, arbitrary, quickCheck, )

import Control.Monad (liftM2, liftM3, guard, )

import qualified Data.Array.Accelerate.Data.Complex as Complex
import Data.Complex (Complex((:+)), cis, )


tolerance :: Double
tolerance = 1e-10

approxEqualAbs ::
   (A.RealFloat a) =>
   Exp a -> Exp a -> Exp a -> Exp Bool
approxEqualAbs eps x y =
   abs (x-y) A.<= eps

approxEqualComplexAbs ::
   (A.RealFloat a) =>
   Exp a -> Exp (Complex a) -> Exp (Complex a) -> Exp Bool
approxEqualComplexAbs eps x y =
   Complex.magnitude (x-y) A.<= eps


genComplex :: QC.Gen (Complex Double)
genComplex = liftM2 (:+) (QC.choose (-1,1)) (QC.choose (-1,1))

newtype Normed0 = Normed0 (Exp (Complex Double))
   deriving (Show)

instance Arbitrary Normed0 where
   arbitrary = fmap (Normed0 . A.constant) genComplex


data Normed1 = Normed1 Int (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1 where
   arbitrary =
      fmap
         (\xs ->
            let len = length xs
            in  Normed1 len $ A.use $ A.fromList (Z :. len) xs) $
      QC.listOf genComplex


floorPowerOfTwo :: Int -> Int
floorPowerOfTwo len =
   2 ^ (floor (logBase 2 (fromIntegral len :: Double)) :: Int)

data Normed1PowerTwo = Normed1PowerTwo Int (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1PowerTwo where
   arbitrary =
      fmap
         (\xs ->
            let len = floorPowerOfTwo $ length xs
            in  Normed1PowerTwo len $ A.use $ A.fromList (Z :. len) xs) $
      liftM2 (:) genComplex (QC.listOf genComplex)


{-
For every reasonable pair of powers of 3 and 5
it computes the largest power of 2,
such that their product is below @n@.
-}
floor5Smooth :: Int -> Int
floor5Smooth n =
   fromInteger $
   maximum $ map (maximum . floorSmooths 2 5 (fromIntegral n)) $
   floorSmooths 2 3 (fromIntegral n) $
   fromIntegral $ floorPowerOfTwo n

{- |
@floorSmooths a b n m@
replaces successively @a@ factors in @m@ by @b@ factors
while keeping the product below @n@.
-}
floorSmooths :: Integer -> Integer -> Integer -> Integer -> [Integer]
floorSmooths a b n =
   let divMany k =
          case divMod k a of
             (q,r) -> guard (r==0) >> if q>n then divMany q else go q
       go m  =  m : divMany (m*b)
   in  go


data Normed1Smooth5 = Normed1Smooth5 Int (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1Smooth5 where
   arbitrary =
      fmap
         (\xs ->
            let len = floor5Smooth $ length xs
            in  Normed1Smooth5 len $ A.use $ A.fromList (Z :. len) xs) $
      liftM2 (:) genComplex (QC.listOf genComplex)


data Normed1Even = Normed1Even Int (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1Even where
   arbitrary =
      fmap
         (\xs ->
            let len = 2 * length xs
            in  Normed1Even len $ A.use $ A.fromList (Z :. len) $
                concatMap (\(x0,x1) -> [x0,x1]) xs) $
      QC.listOf $ liftM2 (,) genComplex genComplex


data
   Normed1Pair =
      Normed1Pair Int
         (Acc (Array DIM1 (Complex Double)))
         (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1Pair where
   arbitrary =
      fmap
         (\xys ->
            let len = length xys
                (xs, ys) = unzip xys
            in  Normed1Pair len
                   (A.use $ A.fromList (Z :. len) xs)
                   (A.use $ A.fromList (Z :. len) ys)) $
      QC.listOf $ liftM2 (,) genComplex genComplex


data
   Normed1Triple =
      Normed1Triple Int
         (Acc (Array DIM1 (Complex Double)))
         (Acc (Array DIM1 (Complex Double)))
         (Acc (Array DIM1 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed1Triple where
   arbitrary =
      fmap
         (\xyzs ->
            let len = length xyzs
                (xs, ys, zs) = unzip3 xyzs
            in  Normed1Triple len
                   (A.use $ A.fromList (Z :. len) xs)
                   (A.use $ A.fromList (Z :. len) ys)
                   (A.use $ A.fromList (Z :. len) zs)) $
      QC.listOf $ liftM3 (,,) genComplex genComplex genComplex


data Normed2 = Normed2 Int Int (Acc (Array DIM2 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed2 where
   arbitrary = do
      xs <- liftM2 (:) QC.arbitrary $ QC.listOf genComplex
      let len = length xs
      height <- QC.choose (1, round (sqrt (fromIntegral len :: Double)))
      let width = div len height
      return $ Normed2 width height $
         A.use $ A.fromList (Z :. height :. width) xs


data Normed3 = Normed3 Int Int Int (Acc (Array DIM3 (Complex Double)))
   deriving (Show)

instance Arbitrary Normed3 where
   arbitrary = do
      xs <- liftM2 (:) QC.arbitrary $ QC.listOf genComplex
      let len = length xs
          lenThd = round $ (fromIntegral len :: Double) ** recip 3
      height <- QC.choose (1, lenThd)
      width <- QC.choose (1, lenThd)
      let depth = div len (width*height)
      return $ Normed3 width height depth $
         A.use $ A.fromList (Z :. depth :. height :. width) xs


-- duplicate of Private.cycleDim3
cycleDim3 :: (A.Elt a) => Transform DIM3 a
cycleDim3 arr =
   A.backpermute
      (Exp.modify (expr:.expr:.expr:.expr)
         (\(sh:.k:.m:.n) -> (sh:.n:.k:.m)) $
       A.shape arr)
      (Exp.modify (expr:.expr:.expr:.expr)
         (\(ix:.n:.k:.m) -> (ix:.k:.m:.n)))
      arr


basisVector :: Int -> Int -> Array DIM1 (Complex Double)
basisVector len freq =
   A.fromList (Z:.len) $
   map (\k -> cis $ 2*pi * fromIntegral (k*freq) / fromIntegral len) $
   iterate (1+) 0


norm2 ::
   Acc (Array DIM1 (Complex Double)) ->
   Acc (A.Scalar Double)
norm2  =
   A.map sqrt . A.sum . A.map (Exp.modify (expr:+expr) $ \(r:+i) -> r*r+i*i)

scalarProduct ::
   Acc (Array DIM1 (Complex Double)) ->
   Acc (Array DIM1 (Complex Double)) ->
   Acc (A.Scalar (Complex Double))
scalarProduct xs ys =
   A.foldAll (+) 0 $ A.zipWith (*) xs (A.map Complex.conjugate ys)

complexFromReal ::
   Acc (Array DIM1 Double) ->
   Acc (Array DIM1 (Complex Double))
complexFromReal = A.map (A.lift . (:+0))

toSelfAdjoint :: Transform DIM1 (Complex Double)
toSelfAdjoint x =
   A.zipWith (+) x $
   A.map Complex.conjugate $ Cyclic.reverse x



infixl 6 .+.

(.+.) ::
   Acc (Array DIM1 (Complex Double)) ->
   Acc (Array DIM1 (Complex Double)) ->
   Acc (Array DIM1 (Complex Double))
(.+.) = A.zipWith (+)


quickCheckWithSign ::
   (Arbitrary normed, Show normed) =>
   (Planned.Sign Double -> normed -> Bool) -> IO ()
quickCheckWithSign = quickCheck


{-
should be replaced by (A.==) in future
-}
class (A.Shape sh, A.Slice sh) => EqShape sh where
   eqShape :: Exp sh -> Exp sh -> Exp Bool

instance EqShape Z where
   eqShape _ _ = A.constant True

instance (EqShape sh, i ~ Int) => EqShape (sh:.i) where
   eqShape =
      Exp.modify2 (expr:.expr) (expr:.expr) $
      \(sh0:.n0) (sh1:.n1) ->
         n0 A.== n1  A.&&  eqShape sh0 sh1


infix 4 =~=

(=~=) ::
   (EqShape sh) =>
   Acc (Array sh (Complex Double)) ->
   Acc (Array sh (Complex Double)) ->
   Acc (A.Scalar Bool)
(=~=) xs ys =
   A.map (eqShape (A.shape xs) (A.shape ys)  A.&&) $ A.and $ A.flatten $
   A.zipWith (approxEqualComplexAbs (A.constant tolerance)) xs ys


run :: Acc (A.Scalar Bool) -> Bool
run = Acc.the . AI.run


tests :: [(String, IO ())]
tests =
   ("fourier generic vs. preprocessed dit2",
      quickCheck $ \sign (Normed1PowerTwo len x) -> run $
             Planned.transform sign len x
             =~=
             Prep.dit2 sign len x) :
   ("fourier generic vs. preprocessed dif2",
      quickCheck $ \sign (Normed1PowerTwo len x) -> run $
             Planned.transform sign len x
             =~=
             Prep.dif2 sign len x) :
   ("fourier generic vs. preprocessed ditSplitRadix",
      quickCheck $ \sign (Normed1PowerTwo len x) -> run $
             Planned.transform sign len x
             =~=
             Prep.ditSplitRadix sign len x) :
   ("fourier generic vs. adhoc dit2",
      quickCheck $ \sign (Normed1PowerTwo len x) -> run $
             Planned.transform sign len x
             =~=
             Adhoc.dit2 (A.constant sign) x) :
   ("fourier generic vs. adhoc ditSplitRadix",
      quickCheck $ \sign (Normed1PowerTwo len x) -> run $
             Planned.transform sign len x
             =~=
             Adhoc.ditSplitRadix (A.constant sign) x) :
   ("fourier generic vs. adhoc dit235",
      quickCheck $ \sign (Normed1Smooth5 len x) -> run $
             Planned.transform sign len x
             =~=
             Adhoc.dit235 (A.constant sign) x) :
   ("fourier adhoc chirp 2 vs. chirp 235",
      quickCheck $ \sign (Normed1 _len x) -> run $
             Adhoc.transformChirp2 (A.constant sign) x
             =~=
             Adhoc.transformChirp235 (A.constant sign) x) :
   ("fourier generic vs. adhoc auto",
      quickCheck $ \sign (Normed1 len x) -> run $
             Planned.transform sign len x
             =~=
             Adhoc.transform (A.constant sign) x) :
   ("fourier generic vs. adhoc chirp 235",
      quickCheck $ \sign (Normed1 len x) -> run $
             Planned.transform sign len x
             =~=
             Adhoc.transformChirp235 (A.constant sign) x) :
   ("fourier generic vs. chirp2",
      quickCheck $ \sign (Normed1 len x) -> run $
             Planned.transform sign len x
             =~=
             Planned.transformChirp2 sign len x) :
   ("fourier generic vs. chirp235",
      quickCheck $ \sign (Normed1 len x) -> run $
             Planned.transform sign len x
             =~=
             Planned.transformChirp235 sign len x) :
   ("homogeneity",
      quickCheck $ \sign (Normed0 x) (Normed1 len y) -> run $
         let transform = Planned.transform sign len
         in  transform (A.map (x*) y)
             =~=
             A.map (x*) (transform y)) :
   ("additivity",
      quickCheck $ \sign (Normed1Pair len x y) -> run $
         let transform = Planned.transform sign len
         in  A.zipWith (+) (transform x) (transform y)
             =~=
             transform (A.zipWith (+) x y)) :
   ("basis vector",
      quickCheck $ \(Normed1 len _x) kp -> run $
         let transform = Planned.transform Planned.inverse len
             k = mod kp len
             unit =
                A.use $ A.fromList (Z:.len) $
                replicate k 0 ++ 1 : repeat 0
         in  transform unit
             =~=
             A.use (basisVector len k)) :
   ("fourier inverse",
      quickCheck $ \(Normed1 len x) -> run $
             x =~=
             (scaleDown $
              Planned.transform Planned.forward len $
              Planned.transform Planned.inverse len x)) :
   ("double fourier = reverse",
      quickCheck $ \sign (Normed1 len x) -> run $
         let transform = Planned.transform sign len
         in  x =~=
             (Cyclic.reverse $
              scaleDown $
              transform $
              transform x)) :
   ("fourier of reverse",
      quickCheck $ \sign (Normed1 len x) -> run $
         let transform = Planned.transform sign len
         in  Cyclic.reverse (transform x) =~=
             transform (Cyclic.reverse x)) :
   ("fourier of conjugate",
      quickCheck $ \sign (Normed1 len x) -> run $
         let transform = Planned.transform sign len
         in  (A.map Complex.conjugate $ transform x)
             =~=
             (transform $
              A.map Complex.conjugate $ Cyclic.reverse x)) :
   ("isometry",
      quickCheck $ \sign (Normed1 len x) -> run $
         let transform = Planned.transform sign len
         in  A.zipWith
                (approxEqualAbs $ A.constant tolerance)
                (norm2 $ transform x)
                (A.map (A.constant (sqrt (fromIntegral len)) *) $ norm2 x)) :
   ("unitarity",
      quickCheck $ \sign (Normed1Pair len x y) -> run $
         let transform = Planned.transform sign len
         in  A.zipWith
                (approxEqualComplexAbs $ A.constant tolerance)
                (scalarProduct (transform x) (transform y))
                (A.map (A.constant (fromIntegral len) *) $
                 scalarProduct x y)) :
   ("convolution commutativity",
      quickCheck $ \(Normed1Pair len x y) -> run $
         let (.*.) = Planned.convolveCyclic len
         in  x .*. y
             =~=
             y .*. x) :
   ("convolution associativity",
      quickCheck $ \(Normed1Triple len x y z) -> run $
         let (.*.) = Planned.convolveCyclic len
         in  (x .*. y) .*. z
             =~=
             x .*. (y .*. z)) :
   ("convolution distributivity",
      quickCheck $ \(Normed1Triple len x y z) -> run $
         let (.*.) = Planned.convolveCyclic len
         in  x .*. (y .+. z)
             =~=
             (x .*. y) .+. (x .*. z)) :
   ("convolution karatsuba rec vs. loop",
      quickCheck $ \(Normed1 len xy) -> run $
         let x = A.map Complex.real xy
             y = A.map Complex.imag xy
         in  complexFromReal (ConvPrep.karatsuba len x y)
             =~=
             complexFromReal (Convolution.karatsuba x y)) :
{-
   No instance for (A.IsNum (Complex Double))
      arising from a use of 'ConvPrep.karatsuba'
-}
   ("convolution karatsuba",
      quickCheck $ \(Normed1Pair len x y) -> run $
         let resultLen = max 0 $ 2*len-1
         in  Convolution.complex Convolution.karatsuba x y
             =~=
             Planned.convolveCyclic resultLen
                (Sliced.pad 0 (A.constant resultLen) x)
                (Sliced.pad 0 (A.constant resultLen) y)) :
   ("convolution cyclic karatsuba",
      quickCheck $ \(Normed1Pair len x y) -> run $
             Convolution.complex
                (Convolution.cyclic Convolution.karatsuba) x y
             =~=
             Planned.convolveCyclic len
                (Sliced.pad 0 (A.constant len) x)
                (Sliced.pad 0 (A.constant len) y)) :
   ("real to spectrum",
      quickCheck $ \(Normed1Even len x) -> run $
         let xr = A.map Complex.real x
         in  FourierReal.toSpectrum
                (Planned.transform Planned.forward (div len 2)) xr
             =~=
             Planned.transform Planned.forward len (complexFromReal xr)) :
   ("real from spectrum",
      quickCheck $ \(Normed1Even len x) -> run $
         let xSelfAdjoint = toSelfAdjoint x
         in  (complexFromReal $
              FourierReal.fromSpectrum
                 (Planned.transform Planned.inverse (div len 2)) xSelfAdjoint)
             =~=
             Planned.transform Planned.inverse len xSelfAdjoint) :
   ("real to and from spectrum",
      quickCheck $ \(Normed1Even len x) -> run $
         let xr = A.map Complex.real x
             len2 = div len 2
         in  (scaleDown $ complexFromReal $
              FourierReal.fromSpectrum (Planned.transform Planned.inverse len2) $
              FourierReal.toSpectrum (Planned.transform Planned.forward len2) xr)
             =~=
             complexFromReal xr) :
   ("real from and to spectrum",
      quickCheck $ \(Normed1Even len x) -> run $
         let len2 = div len 2
             xSelfAdjoint = toSelfAdjoint x
         in  (scaleDown $
              FourierReal.toSpectrum (Planned.transform Planned.forward len2) $
              FourierReal.fromSpectrum (Planned.transform Planned.inverse len2) $
              xSelfAdjoint)
             =~=
             xSelfAdjoint) :
   ("double real to spectrum, even",
      quickCheck $ \(Normed1Even len x) -> run $
         let xr = A.map Complex.real x
             xi = A.map Complex.imag x
             transform = Planned.transform Planned.forward (div len 2)
             (specr,speci) =
                A.unzip $ FourierReal.untangleSpectra $
                Planned.transform Planned.forward len x
         in  A.zipWith (A.&&)
                (FourierReal.toSpectrum transform xr =~= specr)
                (FourierReal.toSpectrum transform xi =~= speci)) :
   ("double real to spectrum, arbitrary",
      quickCheck $ \(Normed1 len x) -> run $
         let xr = complexFromReal $ A.map Complex.real x
             xi = complexFromReal $ A.map Complex.imag x
             transform = Planned.transform Planned.forward len
             (specr,speci) =
                A.unzip $ FourierReal.untangleSpectra $ transform x
         in  A.zipWith (A.&&)
                (transform xr =~= specr)
                (transform xi =~= speci)) :
   ("entangle and untangle spectrum of real data",
      quickCheck $ \(Normed1Pair _len x y) -> run $
         let (xt,yt) =
                A.unzip $
                A.map (A.uncurry FourierReal.untangleCoefficient) $
                A.zipWith FourierReal.entangleCoefficient x y
         in  A.zipWith (A.&&) (x =~= xt) (y =~= yt)) :
   ("double real from spectrum",
      quickCheck $ \(Normed1 len x) -> run $
         let imagUnit = A.constant (0:+1)
             xSelfAdjoint = toSelfAdjoint x
             ySelfAdjoint = toSelfAdjoint $ A.map (imagUnit*) x
             transform = Planned.transform Planned.inverse len
             (xSignal,ySignal) =
                A.unzip $ FourierReal.twoFromSpectrum transform $
                A.zip xSelfAdjoint ySelfAdjoint
         in  A.zipWith (A.&&)
                (transform xSelfAdjoint =~= A.map (A.lift . (:+0)) xSignal)
                (transform ySelfAdjoint =~= A.map (A.lift . (:+0)) ySignal)) :
   ("transform2d vs. transposition, preprocessed",
      quickCheckWithSign $ \sign (Normed2 width height x) -> run $
         let transformH =
                Prep.transform2d
                   (Prep.SubTransformPair
                      (Planned.transform sign width)
                      (Planned.transform sign height))
             transformV =
                Prep.transform2d
                   (Prep.SubTransformPair
                      (Planned.transform sign height)
                      (Planned.transform sign width))
         in  LinAlg.transpose (transformH x)
             =~=
             transformV (LinAlg.transpose x)) :
   ("transform2d vs. transposition, adhoc",
      quickCheckWithSign $ \sign (Normed2 _width _height x) -> run $
         let transform =
                Adhoc.transform2d
                   (Adhoc.SubTransform
                      (Adhoc.transformChirp2 (A.constant sign)))
         in  LinAlg.transpose (transform x)
             =~=
             transform (LinAlg.transpose x)) :
   ("transform3d vs. transposition, preprocessed",
      quickCheckWithSign $ \sign (Normed3 width height depth x) -> run $
         let transformH =
                Prep.transform3d
                   (Prep.SubTransformTriple
                      (Planned.transform sign width)
                      (Planned.transform sign height)
                      (Planned.transform sign depth))
             transformV =
                Prep.transform3d
                   (Prep.SubTransformTriple
                      (Planned.transform sign height)
                      (Planned.transform sign depth)
                      (Planned.transform sign width))
         in  cycleDim3 (transformH x)
             =~=
             transformV (cycleDim3 x)) :
   ("transform2d vs. transposition, adhoc",
      quickCheckWithSign $ \sign (Normed3 _width _height _depth x) -> run $
         let transform =
                Adhoc.transform2d
                   (Adhoc.SubTransform
                      (Adhoc.transformChirp2 (A.constant sign)))
         in  LinAlg.transpose (transform x)
             =~=
             transform (LinAlg.transpose x)) :
   []
