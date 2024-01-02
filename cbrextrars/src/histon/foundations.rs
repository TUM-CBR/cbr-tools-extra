use std::any::{Any, TypeId};
use std::{rc::Rc, collections::LinkedList};

use super::extend::Extend;

enum RelationError {
    IncorrectColumnType { column : String, type_id : TypeId },
    IncorrectColumnCount { expected : usize, actual : usize }
}

type RelationResult<TResult> = Result<TResult, RelationError>;

impl RelationError {

    pub fn raise_incorrect_column_count<T>(
        expected : usize,
        actual : usize
    ) -> RelationResult<T> {

        Result::Err(Self::IncorrectColumnCount{ expected, actual })
    }

    pub fn raise_incorrect_column_type<T>(
        column : String,
    ) -> RelationResult<T>
    where T : Any {
        return Result::Err(RelationError::IncorrectColumnType { column, type_id: TypeId::of::<T>() })
    }

    pub fn incorrect_column_type<T>(
        column : String,
    ) -> RelationError
    where T : Any {
        return RelationError::IncorrectColumnType { column, type_id: TypeId::of::<T>() }
    }
}

struct ToArgsIterator<TValue> {
    arg : TValue
}

impl<TValue> Iterator for ToArgsIterator<TValue> {
    type Item = TValue;

    fn next(&mut self) -> Option<Self::Item> {
        panic!("not implemented")
    }
}

impl<TValue> FromIterator<TValue> for ToArgsIterator<TValue> {

    fn from_iter<T: IntoIterator<Item = TValue>>(iter: T) -> Self {
        panic!("not implemented")
    }
}

pub trait ToArgs : Any + Clone {
    
    fn to_args(
        columns : &Vec<String>,
        args : &Vec<Box<dyn Any>>
    ) -> Result<ToArgsIterator<Self>, RelationError>;
}

pub struct SelectIterator<Values> {
    pub values : Vec<Values>
}

impl<TValue> FromIterator<TValue> for SelectIterator<TValue> {

    fn from_iter<T: IntoIterator<Item = TValue>>(iter: T) -> Self {
        panic!("not implemented")
    }
}

pub trait SelectDispatchFn<Args, TResult> {

    fn dispatch(
        &self,
        columns : &Vec<String>,
        args: &Vec<Box<dyn Any>>
    ) -> RelationResult<SelectIterator<TResult>>;
}

impl<A1> ToArgs for (A1,)
    where A1 : Any + Clone {

    fn to_args(
        columns : &Vec<String>,
        args : &Vec<Box<dyn Any>>) -> RelationResult<ToArgsIterator<(A1,)>> {

        if args.len() != 1 {
            return RelationError::raise_incorrect_column_count(1,columns.len())
        }

        return args[0].downcast::<ToArgsIterator<A1>>()
            .map(|v| {
                v.map(|v| { (v.clone(),) }).collect()
            })
            .or_else(|e| RelationError::raise_incorrect_column_type(columns[0]));
    }
}

impl<A1,A2> ToArgs for (A1,A2)
    where
        A1 : Any + Clone
        , A2 : Any + Clone {

    fn to_args(
        columns : &Vec<String>,
        args : &Vec<Box<dyn Any>>
    ) -> RelationResult<ToArgsIterator<(A1,A2)>> {

        if args.len() != 2 {
            return RelationError::raise_incorrect_column_count(2, columns.len())
        }

        let arg1 = args[0]
            .downcast::<ToArgsIterator<A1>>()
            .map_err(|e| { RelationError::incorrect_column_type::<Box<A1>>(columns[0]) });
        
        let arg2 =
            args[1].downcast::<ToArgsIterator<A2>>()
            .map_err(|b| { RelationError::incorrect_column_type::<Box<A2>>(columns[1])});

        return arg1.and_then(|a1| { 
            arg2.map(|a2|{
                a1.zip(a2).collect()
            })
        });
    }
}

impl<F, Args, TResult> SelectDispatchFn<Args, TResult> for F
    where
        F : FnMut(&Args) -> TResult
        , Args : ToArgs {
    
    fn dispatch(
        &self,
        columns: &Vec<String>,
        args: &Vec<Box<dyn Any>>
    ) -> RelationResult<SelectIterator<TResult>> {
        
        let res = ToArgs::to_args(columns, args).map(|args: Box<Args>|{ self(args.as_ref())});
        
    }
}

pub trait Relation {
    fn select<F, Args, TResult>(
        &self,
        columns : &Vec<String>,
        select : F
    ) -> SelectIterator<TResult>
    where F : SelectDispatchFn<Args, TResult>;
}