{-# LANGUAGE AllowAmbiguousTypes    #-}
module Bio.Cbr.Msa.Algorithms where
    import Conduit ((.|), ConduitT, foldlC, runConduitRes, mapAccumWhileC, sinkNull, sourceHandle)
    import Control.Lens (ASetter', Lens, Lens', view)
    import Control.Lens.Combinators (assign, modifying, over, use)
    import Control.Monad.Reader.Class (MonadReader)
    import Control.Monad.State (evalState, get, put)
    import Control.Monad.State.Strict (execStateT)
    import Control.Monad.State.Class (MonadState)
    import qualified Data.Conduit.Combinators as Conduit
    import qualified Data.Conduit.Text as TextC
    import Data.Foldable (foldr)
    import Data.Generics.Labels (Field, Field', fieldLens)
    import Data.Map (Map)
    import qualified Data.Map as Map
    import Data.Maybe (fromMaybe)
    import Data.Proxy (Proxy)
    import Data.Text (Text)
    import qualified Data.Text as Text
    import Data.Vector (Vector, snoc)
    import qualified Data.Vector as Vector
    import GHC.Generics (Generic)
    import GHC.OverloadedLabels (IsLabel)
    import GHC.TypeLits (Symbol)
    import Prelude hiding (map, foldl, foldr)
    import Control.Lens.Internal.Context (IndexedComonadStore(context))
    import qualified Data.Conduit.Text as Text
    import qualified Control.Monad.State as Control.Monad.Trans
    
    data Indexed a = Indexed { unIndex :: Int, unValue :: a}

    type IText = Indexed Text

    type IGapsCount = Map Int (Vector Int)

    type IGapsByPosition = Map Int Int

    type IGapNegligenceScore = Map Int Int

    data PipelineContext state computed =
        PipelineContext
        { unState :: state
        , unComputed :: state -> computed
        } deriving Generic

    class Empty value where
        empty :: value

    instance Empty (Map k v) where
        empty = Map.empty

    instance Empty a => Empty (value -> a) where
        empty = const empty

    instance Empty (Vector a) where
        empty = Vector.empty

    instance (Empty state, Empty synth) => Empty (PipelineContext state synth) where
        empty = PipelineContext empty empty

    iAssign :: (MonadState (PipelineContext state synth) m) => Lens' state (Map Int value) -> Indexed a -> value -> m ()
    iAssign theLens key value = modifying (fieldLens @"unState" . theLens) (Map.insert (unIndex key) value)

    iSynthAssign :: (MonadState (PipelineContext state synth) m) => Lens' synth (Map Int value) -> Indexed a -> (state -> value) -> m ()
    iSynthAssign theLens key value = modifying (fieldLens @"unComputed") updater
        where
            updater synth context = over theLens (Map.insert (unIndex key) (value context)) (synth context)

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
        MsaField "lGapsCount" s = Field' "lGapsCount" s IGapsCount
        MsaField "lGapsNegligence" s = Field' "lGapsNegligence" s IGapNegligenceScore

    mapAccum :: Monad m => (t -> a -> (a, b)) -> a -> ConduitT t b m a
    mapAccum mapper = mapAccumWhileC accumulator
        where
            accumulator input = Right . mapper input

    indexedC :: Monad m => ConduitT a (Indexed a) m ()
    indexedC = () <$ mapAccum withIndex 0
        where
            withIndex value state = (state + 1, Indexed { unIndex = state, unValue = value})

    identifyGapsM :: 
        (MonadState (PipelineContext state computed) m, MsaField "lGapsCount" state) =>
        IText ->
        m IText
    identifyGapsM input = do
        let gaps = identifyGaps $ unValue input
        iAssign (fieldLens @"lGapsCount") input gaps
        pure input

    identifyGapsC
        :: (MonadState (PipelineContext state computed) m, MsaField "lGapsCount" state)
        => ConduitT IText IText m ()
    identifyGapsC = Conduit.mapM identifyGapsM

    countByPosition :: IGapsCount -> IGapsByPosition
    countByPosition = Map.foldr acc Map.empty
        where
            acc gaps result = foldr (Map.insertWith (+) 1) result gaps

    scoreByGapNegligenceC ::
        forall state synth m .
        ( MonadState (PipelineContext state synth) m
        , MsaField "lGapsNegligence" synth
        , MsaField "lGapsCount" state
        ) =>
        ConduitT IText IText m ()
    scoreByGapNegligenceC = Conduit.mapM $ \iText -> do
        iSynthAssign (fieldLens @"lGapsNegligence") iText $ foldTextWithIndex agg (const 0) (unValue iText)
        pure iText

        where
            gapsLens :: state -> IGapsByPosition
            gapsLens = countByPosition . view (fieldLens @"lGapsCount")
            agg :: (state -> Int) -> Int -> Char -> state -> Int
            agg s i c
                -- If there is a gap at a position, there must be an entry
                -- in the map. Otherwise there is a bug
                | isGap c = \state -> gapsLens state Map.! i + s state
                | otherwise = s

    newtype ScoreByGapNegligenceState =
        ScoreByGapNegligenceState
        {
            lGapsCount :: IGapsCount
        } deriving Generic

    instance Empty ScoreByGapNegligenceState where
        empty = ScoreByGapNegligenceState empty

    newtype ScoreByGapNegligenceSynth =
        ScoreByGapNegligenceSynth
        {
            lGapsNegligence :: IGapNegligenceScore
        } deriving Generic

    instance Empty ScoreByGapNegligenceSynth where
        empty = ScoreByGapNegligenceSynth empty

    scoreByGapNegigencePipeline source =
        sourceHandle source
        .| TextC.decode TextC.utf8
        .| TextC.lines
        .| indexedC
        .| identifyGapsC
        .| scoreByGapNegligenceC
        .| sinkNull

    scoreByGapNegigence source = do
        result <- flip execStateT empty $ runConduitRes $ scoreByGapNegigencePipeline source
        pure $ fixT result

        where
            fixT :: PipelineContext ScoreByGapNegligenceState ScoreByGapNegligenceSynth -> ()
            fixT _ = ()
