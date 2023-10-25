module Bio.Cbr.IO where

    import Control.Concurrent.MVar (MVar, modifyMVar, newMVar, takeMVar)
    import Control.Concurrent.Chan (Chan, newChan)
    import Control.Monad (foldM)
    import Data.ByteString (ByteString)
    import qualified Data.ByteString.Lazy as ByteString
    import Data.ByteString.Builder as Builder
    import Data.Map (Map)
    import Data.Word (Word8)
    import qualified Data.Map as Map
    import Foreign.Ptr (Ptr, plusPtr)
    import Foreign.Storable (peekElemOff)
    import System.IO.Unsafe (unsafePerformIO)
    
    type PyHandle = Int

    newtype PyStream = PyStream { channel :: Chan ByteString }

    type Handles = MVar (Map Int PyStream)

    {-# NOINLINE handles #-}
    handles :: Handles
    handles = unsafePerformIO $ newMVar Map.empty

    foreign export ccall "newPyStream" newPyStream :: IO Int

    newPyStream :: IO PyHandle
    newPyStream = modifyMVar handles $ \values -> do
        let key = (+1) . maybe 0 fst $ Map.lookupMax values
        channel <- newChan
        pure (Map.insert key PyStream{channel = channel} values, key)

    fromPointer :: Ptr Word8 -> Int -> IO ByteString
    fromPointer ptr count =
        ByteString.toStrict . Builder.toLazyByteString <$> foldM acc (Builder.byteString "") [0 .. count]
        where
            acc builder offset = (builder <>) . word8 <$> peekElemOff ptr offset

    foreign export ccall "writePyStream" writePyStream :: Ptr Word8 -> Int -> PyHandle -> IO ()

    writePyStream :: Ptr Word8 -> Int -> PyHandle -> IO ()
    writePyStream bytes count handle = do
        bs <- fromPointer bytes count
        stream <- (Map.! handle) <$> takeMVar handles
        pure ()