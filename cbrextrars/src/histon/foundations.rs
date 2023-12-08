use core::iter::Iterator;
use core::slice::Iter;
use std::any::Any;

pub trait IsTuple {
    fn get<T : Any>(&self, column: &String) -> Option<&T>;
}

pub trait IndexInstance<Key> where Key : Eq {
    type Iterator : Iterator<Item = u64>;

    // Get a slice of values ordered as implied by the index. This has the format
    // from, to are inclusive meaning that from == to gets all the values deemed
    // equal according to the index. Ommiting either returns all other elements.
    // This implies that from(None, None) should return all indices.
    fn slice(&self, from: Option<&Key>, to: Option<&Key>) -> Self::Iterator;
}

pub struct IndexIterator<'a, Tuple> where Tuple : IsTuple {
    items : &'a mut Iter<'a, RelationRow<Tuple>>
}

impl<'a, Tuple> Iterator for IndexIterator<'a, Tuple>
    where
        Tuple : IsTuple {

    type Item = &'a RelationRow<Tuple>;

    fn next(&mut self) -> Option<Self::Item> {
        return self.items.next();
    }
}

pub trait Index<Key> where Key : Eq {
    type Instance : IndexInstance<Key>;

    fn create_index<Tuple>(items: &mut IndexIterator<Tuple>) -> Self::Instance
        where Tuple : IsTuple;
}

pub struct RelationRow<Tuple> {
    pub row : u64,
    pub tuple : Tuple
}

pub trait Relation<'t> {
    type Tuple : IsTuple + 't;
    type Iterator : Iterator<Item = &'t RelationRow<Self::Tuple>>;

    fn tuples(&'t self) -> Self::Iterator;
}

pub trait IndexedRelation<'t, Key> where Key : Eq {
    type Index : IndexInstance<Key>;
    type Out : Relation<'static>;

    fn slice(&self, from : Option<&Key>, to : Option<&Key>) -> Self::Out;
}

pub trait ToIndexed {
    type Out<Key> : IndexedRelation<'static, Key> where Key : Eq;

    fn indexed<I, Key>(&self, index: &I) -> Self::Out<Key>
        where
            I : Index<Key>
            , Key : Eq;
}