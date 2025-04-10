use nom::{
    number::complete::{le_f64, le_i32},
    IResult,
};

#[derive(Debug, PartialEq)]
pub struct Summary {
    pub start_epoch: f64,
    pub end_epoch: f64,
    pub target: i32,
    pub center: i32,
    pub frame_id: i32,
    pub data_type: i32,
    pub initial_addr: i32,
    pub final_addr: i32,
}

impl Summary {
    pub fn parse(input: &[u8]) -> IResult<&[u8], Self> {
        let (input, start_epoch) = le_f64(input)?;
        let (input, end_epoch) = le_f64(input)?;

        let (input, target) = le_i32(input)?;
        let (input, center) = le_i32(input)?;
        let (input, frame_id) = le_i32(input)?;
        let (input, data_type) = le_i32(input)?;
        let (input, initial_addr) = le_i32(input)?;
        let (input, final_addr) = le_i32(input)?;
        Ok((
            input,
            Summary {
                start_epoch,
                end_epoch,
                target: target,
                center: center,
                frame_id: frame_id,
                data_type: data_type,
                initial_addr: initial_addr,
                final_addr: final_addr,
            },
        ))
    }
}
