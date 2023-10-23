module Bio.Cbr where

    foreign export ccall "demo" demo :: IO ()

    demo = print "hello"