-- cf. synthesizer-core:Synthesizer.Generic.Cyclic
{- |
Some small size convolutions using the Karatsuba trick.
We do not use Toom-3 multiplication,
because this requires division by 2 and 6,
and thus 'Fractional' constraints.
-}
module Data.Array.Accelerate.Convolution.Small where


type Pair a = (a,a)

convolvePair ::
   (Num a) =>
   Pair a -> Pair a -> Pair a
convolvePair a b =
   snd $ sumAndConvolvePair a b

sumAndConvolvePair ::
   (Num a) =>
   Pair a -> Pair a -> ((a,a), Pair a)
sumAndConvolvePair (a0,a1) (b0,b1) =
   let sa01 = a0+a1
       sb01 = b0+b1
       ab0ab1 = a0*b0+a1*b1
   in  ((sa01, sb01), (ab0ab1, sa01*sb01-ab0ab1))


type Triple a = (a,a,a)

convolveTriple ::
   (Num a) =>
   Triple a -> Triple a -> Triple a
convolveTriple a b =
   snd $ sumAndConvolveTriple a b

sumAndConvolveTriple ::
   (Num a) =>
   Triple a -> Triple a -> ((a,a), Triple a)
sumAndConvolveTriple (a0,a1,a2) (b0,b1,b2) =
   let ab0 = a0*b0
       dab12 = a1*b1 - a2*b2
       sa01 = a0+a1; sb01 = b0+b1; tab01 = sa01*sb01 - ab0
       sa02 = a0+a2; sb02 = b0+b2; tab02 = sa02*sb02 - ab0
       sa012 = sa01+a2
       sb012 = sb01+b2

       d0 = sa012*sb012 - tab01 - tab02
       d1 = tab01 - dab12
       d2 = tab02 + dab12
   in  ((sa012, sb012), (d0, d1, d2))


type Quadruple a = (a,a,a,a)

convolveQuadruple ::
   (Num a) =>
   Quadruple a -> Quadruple a -> Quadruple a
convolveQuadruple a b =
   snd $ sumAndConvolveQuadruple a b

sumAndConvolveQuadruple ::
   (Num a) =>
   Quadruple a -> Quadruple a -> ((a,a), Quadruple a)
sumAndConvolveQuadruple (a0,a1,a2,a3) (b0,b1,b2,b3) =
   let ab0 = a0*b0
       ab1 = a1*b1
       sa01 = a0+a1; sb01 = b0+b1
       ab01 = sa01*sb01 - (ab0+ab1)
       ab2 = a2*b2
       ab3 = a3*b3
       sa23 = a2+a3; sb23 = b2+b3
       ab23 = sa23*sb23 - (ab2+ab3)
       c0 = ab0  + ab2 - (ab1 + ab3)
       c1 = ab01 + ab23
       ab02 = (a0+a2)*(b0+b2)
       ab13 = (a1+a3)*(b1+b3)
       sa0123 = sa01+sa23
       sb0123 = sb01+sb23
       ab0123 = sa0123*sb0123 - (ab02+ab13)
       d0 = ab13   + c0
       d1 =          c1
       d2 = ab02   - c0
       d3 = ab0123 - c1
   in  ((sa0123, sb0123), (d0, d1, d2, d3))
