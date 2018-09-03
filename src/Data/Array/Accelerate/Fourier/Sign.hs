{-# LANGUAGE ConstraintKinds       #-}
{-# LANGUAGE DeriveDataTypeable    #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies          #-}
module Data.Array.Accelerate.Fourier.Sign where

import           Data.Array.Accelerate.Data.Complex (Complex ((:+)))

import           Data.Array.Accelerate              (Lift (lift), Plain,
                                                     Unlift (unlift))
import qualified Data.Array.Accelerate              as A
import           Data.Array.Accelerate.Array.Sugar  (Elt (eltType, fromElt, toElt),
                                                     EltRepr,
                                                     Tuple (NilTup, SnocTup))
import           Data.Array.Accelerate.Product      (IsProduct (ProdRepr, fromProd, prod, toProd),
                                                     ProdR (ProdRsnoc, ProdRunit),
                                                     TupleIdx (ZeroTupIdx))
import           Data.Array.Accelerate.Smart        (Exp (Exp),
                                                     PreExp (Prj, Tuple))

import           Data.Typeable                      (Typeable)

import qualified Test.QuickCheck                    as QC


newtype Sign a = Sign {getSign :: a}
   deriving (Eq, Show, Typeable)

type instance EltRepr (Sign a) = EltRepr a

instance Elt a => Elt (Sign a) where
   eltType = eltType . getSign
   toElt   = Sign . toElt
   fromElt = fromElt . getSign

instance (cst a) => IsProduct cst (Sign a) where
   type ProdRepr (Sign a) = ((), a)
   fromProd _ (Sign a) = ((), a)
   toProd _ ((), a) = Sign a
   prod _ (Sign _) = ProdRsnoc ProdRunit

instance (Lift Exp a, Elt (Plain a)) => Lift Exp (Sign a) where
   type Plain (Sign a) = Sign (Plain a)
   lift (Sign a) = Exp $ Tuple (NilTup `SnocTup` lift a)

instance Elt a => Unlift Exp (Sign (Exp a)) where
   unlift e = Sign $ Exp $ ZeroTupIdx `Prj` e


forward, inverse :: Num a => Sign a
forward = Sign (-1)
inverse = Sign 1

forwardExp, inverseExp :: (A.Num a) => Exp (Sign a)
forwardExp = lift $ Sign (-1 :: (A.Num a) => Exp a)
inverseExp = lift $ Sign ( 1 :: (A.Num a) => Exp a)

toSign :: (Elt a) => Exp (Sign a) -> Exp a
toSign = getSign . unlift

cis ::
   (A.RealFloat a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp a -> Exp (Complex a)
cis sign w  =  A.lift $ cos w :+ toSign sign * sin w

cisRat ::
   (A.RealFloat a, A.FromIntegral Int a, A.Elt (Complex a)) =>
   Exp (Sign a) -> Exp Int -> Exp Int -> Exp (Complex a)
cisRat sign denom numer =
   cis sign $ 2*pi * A.fromIntegral numer / A.fromIntegral denom


instance (Num a) => QC.Arbitrary (Sign a) where
   arbitrary = QC.elements [forward, inverse]
