use core::array::*;
use std::{any::Any, collections::HashMap};
use super::foundations::*;

struct StaticRelation {
    columns : HashMap<String, Box<dyn Any>>
}

impl Relation for StaticRelation {

    fn select<FOut, Args, TResult>(
        &self,
        columns : &Vec<String>,
        select : FOut
    ) -> SelectIterator<TResult>
    where
        FOut : SelectDispatchFn<Args, TResult> {

    
    let args =
        columns.iter()
        .map(|col| { self.columns[col][0] })
        .collect();
    
    select.dispatch(columns, &args);

    return SelectIterator {
        values : [].to_vec()
    }
}
}