{-# LANGUAGE AllowAmbiguousTypes, ScopedTypeVariables #-}
module Bio.Cbr.Msa.Algorithms where
    import Conduit (ConduitT)
    import Control.Lens (ASetter', Lens')
    import Control.Lens.Combinators (assign, modifying)
    import Control.Monad.State (evalState, get, put)
    import Control.Monad.State.Class (MonadState)
    import qualified Data.Conduit.Combinators as Conduit
    import Data.Foldable (foldr)
    import Data.Generics.Product.Fields (HasField'(..))
    import Data.Map (Map)
    import qualified Data.Map as Map
    import Data.Proxy (Proxy)
    import Data.Text (Text)
    import qualified Data.Text as Text
    import Data.Vector (Vector, empty, snoc)
    import GHC.OverloadedLabels (IsLabel)
    import Prelude hiding (map, foldl, foldr)

    data Indexed a = Indexed { unIndex :: Int, unValue :: a}

    type IText = Indexed Text

    type IGapsCount = Map Int (Vector Int)

    type IGapsByPosition = Map Int Int

    iAssign :: MonadState s m => ASetter' s (Map Int value) -> Indexed a -> value -> m ()
    iAssign theLens key value = modifying theLens (Map.insert (unIndex key) value)


    textFoldlM :: Monad m => (a -> Char -> m a) -> a -> Text -> m a
    textFoldlM acc initial = Text.foldl acc' (pure initial)
        where
            acc' s v = s >>= flip acc v

    identifyGaps :: Text -> Vector Int
    identifyGaps = flip evalState (0 :: Int) . textFoldlM acc empty
        where
            acc s c = do
                i <- get
                put (i + 1)
                pure $ if
                    | c `elem` ['-'] -> snoc s i
                    | otherwise -> s

    identifyGapsM ::  forall s m . (MonadState s m, HasField' "lGapsCount" s IGapsCount) =>  IText -> m IText
    identifyGapsM input = do
        let gaps = identifyGaps $ unValue input
        iAssign (field' @"lGapsCount") input gaps
        pure input

    identifyGapsC :: forall s m . (MonadState s m, HasField' "lGapsCount" s IGapsCount) => ConduitT IText IText m ()
    identifyGapsC = Conduit.mapM identifyGapsM

    countByPosition :: IGapsCount -> IGapsByPosition
    countByPosition = Map.foldr acc Map.empty
        where
            acc gaps result = foldr (Map.insertWith (+) 1) result gaps

