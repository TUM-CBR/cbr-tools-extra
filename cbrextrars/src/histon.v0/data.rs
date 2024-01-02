use std::any::Any;
use std::collections::HashMap;
use std::slice::Iter;

use super::foundations::*;

pub struct HashTuple {
    items : HashMap<String, Box<dyn Any>>
}

impl HashTuple {
    pub fn from_any<T : Any>(column : &String, value: T) -> HashTuple {

        let value: Box<dyn Any> = Box::new(value);
        HashTuple { items: HashMap::from([(column.clone(), value)]) }
    }
}

impl IsTuple for HashTuple {

    fn get<T : Any>(&self, column: &String) -> Option<&T> {
        return self.items[column].downcast_ref::<T>();
    }
}

pub struct VecRelation {
    items : Vec<RelationRow<HashTuple>>
}

impl VecRelation {

    pub fn from_array<T>(
        column: &String,
        values: &Vec<T>
    ) -> VecRelation
    where T : Any + Copy {

        let items = values.into_iter()
            .map(
                |item| { HashTuple::from_any(column, item.clone() )}
            )
            .scan(
                0,
                |state, item| {
                    *state += 1;
                    Some(RelationRow{ row : state.clone() - 1, tuple: item })
                }
            );
        return VecRelation { items: Vec::from_iter(items) }
    }
}

impl<T> From<(String, &Vec<T>)> for VecRelation where T : Any + Copy {
    fn from((column, values) : (String, &Vec<T>)) -> VecRelation {
        return VecRelation::from_array(&column, values)
    }
}

impl<T> From<(&String, &Vec<T>)> for VecRelation where T : Any + Copy {
    fn from((column, values) : (&String, &Vec<T>)) -> VecRelation {
        return (column, values).into()
    }
}

impl<T> From<(&str, &Vec<T>)> for VecRelation where T : Any + Copy {
    fn from((column, values) : (&str, &Vec<T>)) -> VecRelation {
        return (column, values).into()
    }
}

impl<'t> Relation<'t> for VecRelation {

    type Tuple = HashTuple;
    type Iterator = Iter<'t, RelationRow<Self::Tuple>>;

    fn tuples(&'t self) -> Self::Iterator{
        return self.items.iter()
    }
}