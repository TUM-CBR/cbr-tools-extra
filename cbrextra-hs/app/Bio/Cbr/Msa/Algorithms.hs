module Bio.Cbr.Msa.Algorithms where
    import Conduit (ConduitT)
    import Control.Lens (ASetter', Lens')
    import Control.Lens.Combinators (assign, modifying, use)
    import Control.Monad.State (evalState, get, put)
    import Control.Monad.State.Class (MonadState)
    import qualified Data.Conduit.Combinators as Conduit
    import Data.Foldable (foldr)
    import Data.Generics.Product.Fields (HasField'(..))
    import Data.Map (Map)
    import qualified Data.Map as Map
    import Data.Maybe (fromMaybe)
    import Data.Proxy (Proxy)
    import Data.Text (Text)
    import qualified Data.Text as Text
    import Data.Vector (Vector, empty, snoc)
    import GHC.OverloadedLabels (IsLabel)
    import GHC.TypeLits (Symbol)
    import Prelude hiding (map, foldl, foldr)
    
    data Indexed a = Indexed { unIndex :: Int, unValue :: a}

    type IText = Indexed Text

    type IGapsCount = Map Int (Vector Int)

    type IGapsByPosition = Map Int Int

    type IGapNegligenceScore = Map Int (IGapsByPosition -> Int)

    iAssign :: MonadState s m => ASetter' s (Map Int value) -> Indexed a -> value -> m ()
    iAssign theLens key value = modifying theLens (Map.insert (unIndex key) value)

    foldTextLM :: Monad m => (a -> Char -> m a) -> a -> Text -> m a
    foldTextLM acc initial = Text.foldl acc' (pure initial)
        where
            acc' s v = s >>= flip acc v

    forTextM :: Monad m => (Char -> m ()) -> Text -> m ()
    forTextM action = foldTextLM (const action) ()

    foldTextWithIndex :: (a -> Int -> Char -> a) -> a -> Text -> a
    foldTextWithIndex acc initial = flip evalState 0 . foldTextLM acc' initial
        where
            acc' s c = do
                i <- get
                put $ i + 1
                pure $ acc s i c
    
    isGap :: Char -> Bool
    isGap c = c `elem` gapChars
        where
            gapChars = ['-']

    identifyGaps :: Text -> Vector Int
    identifyGaps = foldTextWithIndex acc empty
        where
            acc s i c
                | isGap c = snoc s i
                | otherwise = s

    type family MsaField (name :: Symbol) s where
        MsaField "lGapsCount" s = HasField' "lGapsCount" s IGapsCount
        MsaField "lGapsNegligence" s = HasField' "lGapsNegligence" s IGapNegligenceScore

    identifyGapsM ::  forall s m . (MonadState s m, MsaField "lGapsCount" s) =>  IText -> m IText
    identifyGapsM input = do
        let gaps = identifyGaps $ unValue input
        iAssign (field' @"lGapsCount") input gaps
        pure input

    identifyGapsC :: (MonadState s m, MsaField "lGapsCount" s) => ConduitT IText IText m ()
    identifyGapsC = Conduit.mapM identifyGapsM

    countByPosition :: IGapsCount -> IGapsByPosition
    countByPosition = Map.foldr acc Map.empty
        where
            acc gaps result = foldr (Map.insertWith (+) 1) result gaps

    scoreByGapNegligence ::
        (MonadState s m, MsaField "lGapsNegligence" s) =>
        ConduitT IText IText m ()
    scoreByGapNegligence = Conduit.mapM mapping
        where
            agg :: (IGapsByPosition -> Int) -> Int -> Char -> IGapsByPosition -> Int
            agg s i c
                -- If there is a gap at a position, there must be an entry
                -- in the map. Otherwise there is a bug
                | isGap c = \gapsByPosition -> s gapsByPosition + gapsByPosition Map.! i
                | otherwise = s
            mapping iText = do
                iAssign (field' @"lGapsNegligence") iText $ foldTextWithIndex agg (const 0) (unValue iText)
                pure iText