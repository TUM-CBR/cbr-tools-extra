module Bio.Cbr where

    import Data.Hashable
    foreign export ccall "demo" demo :: IO ()

    demo = print "hello"