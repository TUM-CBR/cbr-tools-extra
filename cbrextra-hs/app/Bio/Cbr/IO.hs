module Bio.Cbr.IO where

    import Control.Concurrent.MVar (MVar, modifyMVar, newMVar, takeMVar, withMVar)
    import Control.Concurrent.Chan (Chan, newChan, readChan, writeChan)
    import Control.Monad (foldM)
    import Data.ByteString (ByteString, packCString)
    import Data.ByteString.Char8 (unpack)
    import Data.ByteString.Builder as Builder
    import Data.Map (Map)
    import Data.Word (Word8)
    import qualified Data.Map as Map
    import Foreign.C.String (CString, newCString)
    import Foreign.Ptr (Ptr)
    import Foreign.Storable (peekElemOff)
    import System.IO.Unsafe (unsafePerformIO)
    
    newtype PyHandle = PyHandle {handleId :: Int}

    newtype PyStream = PyStream { channel :: Chan ByteString }

    type Handles = MVar (Map Int PyStream)

    {-# NOINLINE handles #-}
    handles :: Handles
    handles = unsafePerformIO $ newMVar Map.empty

    foreign export ccall "newPyStream" newPyStream :: IO PyHandle

    newPyStream :: IO PyHandle
    newPyStream = modifyMVar handles $ \values -> do
        let key = (+1) . maybe 0 fst $ Map.lookupMax values
        channel <- newChan
        pure (Map.insert key PyStream{channel = channel} values, PyHandle key)

    foreign export ccall "writePyStream" writePyStream :: CString -> PyHandle -> IO ()

    writePyStream :: CString -> PyHandle -> IO ()
    writePyStream bytes PyHandle{..} = do
        bs <- packCString bytes
        PyStream{..} <- (Map.! handleId) <$> takeMVar handles
        writeChan channel bs

    foreign export ccall "readPyStream" readPyStream :: PyHandle -> IO CString

    readPyStream :: PyHandle -> IO CString
    readPyStream PyHandle{..} = do
        PyStream{..} <- withMVar handles $ pure . (Map.! handleId)
        readChan channel >>= newCString . unpack