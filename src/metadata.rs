use std::env;

use chrono::prelude::*;
use indexmap::IndexMap;
use sysinfo::{Pid, System};
use users::{get_current_gid, get_current_uid, get_group_by_gid, get_user_by_uid};

pub fn collect_starting_metadata(md: &mut IndexMap<&'static str, String>) {
    // Get command line
    let command_line = env::args_os().fold(String::new(), |mut s, a| {
        let arg = a.to_string_lossy();
        if !s.is_empty() {
            s.push(' ')
        };
        s.push_str(&arg);
        s
    });
    md.insert("CommandLine", command_line);

    // Get system info
    let mut sys = System::new();
    sys.refresh_memory();
    sys.refresh_processes();
    if let Some(s) = System::host_name() {
        md.insert("node_name", s);
    }
    md.insert("total_memory", format!("{} K", sys.total_memory()));
    md.insert("physical_cores", format!("{}", num_cpus::get_physical()));
    md.insert("virtual_cores", format!("{}", num_cpus::get()));
    if let Some(proc) = sys.process(Pid::from_u32(unsafe { libc::getpid() as u32 })) {
        md.insert("pid", format!("{}", proc.pid()));
        md.insert("cwd", format!("{:?}", proc.cwd()));
    }
    // Get user and group
    let uid = get_current_uid();
    let gid = get_current_gid();
    if let (Some(user), Some(group)) = (get_user_by_uid(uid), get_group_by_gid(gid)) {
        md.insert(
            "user_group",
            format!(
                "{}:{} ({}:{})",
                user.name().to_string_lossy(),
                group.name().to_string_lossy(),
                uid,
                gid
            ),
        );
    }

    // Get starting time
    let local = Local::now();
    md.insert("starting_data_time", local.to_rfc2822());
}
