module Main where

import qualified Test.Data.Array.Accelerate.Fourier as Fourier

import Data.Tuple.HT (mapFst, )


prefix :: String -> [(String, IO ())] -> [(String, IO ())]
prefix msg =
   map (mapFst (\str -> msg ++ "." ++ str))

main :: IO ()
main =
   mapM_ (\(msg,io) -> putStr (msg++": ") >> io) $
   prefix "Fourier" Fourier.tests
