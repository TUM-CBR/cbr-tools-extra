use super::foundations::*;

pub struct Extend<TRelIn, F> {
    relation : TRelIn,
    columns : Vec<String>,
    f : F
}

impl<TRelIn, F> Relation for Extend<TRelIn, F> where TRelIn : Relation {

    fn select<FOut, Args, TResult>(
            columns : Vec<String>,
            select : FOut
        ) -> SelectIterator<TResult>
        where FOut : FnMut(Args) -> TResult {

        

        return SelectIterator {
            values : [].to_vec()
        }
    }
}
