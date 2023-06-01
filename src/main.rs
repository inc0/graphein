use clap;
use pdbtbx::*;
use std::collections::HashMap;
use petgraph::{graph::Graph, graph::NodeIndex};
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::prelude::*;
use anyhow::{Result, bail};
use log::{debug, info, warn};
use glob::glob;
use rayon::prelude::*;


fn van_der_waals_radius(element: &Element) -> f64 {
    match element {
        Element::C => 1.70,
        Element::H => 1.20,
        Element::N => 1.55,
        Element::O => 1.52,
        Element::P => 1.80,
        Element::S => 1.80,
        Element::Ca => 2.31,
        Element::K => 2.75,
        Element::Na => 2.27,
        Element::Cl => 1.75,
        Element::Mg => 1.73,
        Element::I => 1.98,
        Element::Se => 1.90,
        Element::Cu => 1.40,
        Element::F => 1.47,
        Element::Br => 1.85,
        _ => 1.8, 
    }
}

fn atomic_number(element: &Element) -> u8 {
    match element {
        Element::H => 1,
        Element::C => 6,
        Element::N => 7,
        Element::O => 8,
        Element::F => 9,
        Element::Na => 11,
        Element::Mg => 12,
        Element::P => 15,
        Element::S => 16,
        Element::Cl => 17,
        Element::K => 19,
        Element::Ca => 20,
        Element::Cu => 29,
        Element::Br => 35,
        Element::Se => 34,
        Element::I => 53,
        _ => 0, // Let's assume that 0 means 'unknown'
    }
}

fn valence_electrons(element: &Element) -> u8 {
    match element {
        Element::H => 1,
        Element::C => 4,
        Element::N => 5,
        Element::O => 6,
        Element::F => 7,
        Element::Na => 1,
        Element::Mg => 2,
        Element::P => 5,
        Element::S => 6,
        Element::Cl => 7,
        Element::K => 1,
        Element::Ca => 2,
        Element::Cu => 1, // Typically in a +2 oxidation state, it loses 2 electrons
        Element::Br => 7,
        Element::Se => 6,
        Element::I => 7,
        _ => 0, // Let's assume that 0 means 'unknown'
    }
}

fn electronegativity(element: &Element) -> f64 {
    match element {
        Element::H => 2.20,
        Element::C => 2.55,
        Element::N => 3.04,
        Element::O => 3.44,
        Element::F => 3.98,
        Element::Na => 0.93,
        Element::Mg => 1.31,
        Element::P => 2.19,
        Element::S => 2.58,
        Element::Cl => 3.16,
        Element::K => 0.82,
        Element::Ca => 1.00,
        Element::Cu => 1.90,
        Element::Br => 2.96,
        Element::Se => 2.55,
        Element::I => 2.66,
        _ => 0.0, // Let's assume that 0 means 'unknown'
    }
}


#[derive(Serialize, Deserialize, Debug)]
struct AtomNode {
    id: usize,
    atom_number: u8,
    valence: u8,
    electronegativity: f64,
    charge: isize,
}


fn process_pdb_file(fname: &str, edge_max_dist: &f64) -> Result<()> {
    let (pdb, _errors) = match pdbtbx::open(
        fname,
        StrictnessLevel::Medium
    ) {
        Ok(pdb) => pdb,
        Err(e) => bail!("Error parsing pdb file {} - {:?}", fname, e)
    };

    let tree = pdb.create_atom_rtree();
    let mut protein_graph = Graph::<AtomNode, f64>::new();
    let mut atom_sn_node_id: HashMap<usize, NodeIndex> = HashMap::new();
    
    for atom in pdb.atoms() {
        let ele = match atom.element() {
            Some(e) => e,
            None => continue
        };
        let an = AtomNode {
            id: atom.serial_number(),
            atom_number: atomic_number(ele),
            valence: valence_electrons(ele),
            electronegativity: electronegativity(ele),
            charge: atom.charge(),
        };
        let node_id = protein_graph.add_node(an);
        atom_sn_node_id.insert(atom.serial_number(), node_id);

    }

    for atom in pdb.atoms() {
        let atom_node_id = match atom_sn_node_id.get(&atom.serial_number()) {
            Some(an) => an,
            None => continue
        };
        for neighbor_atom in tree.locate_within_distance(atom.pos(), edge_max_dist * edge_max_dist) {
            let neigh_sn = neighbor_atom.serial_number();
            if atom.pos() == neighbor_atom.pos() {  // Same atom
                continue;
            };
            let node_id =  match atom_sn_node_id.get(&neigh_sn){
                Some(ni) => ni,
                None => continue
            };
            protein_graph.update_edge(*atom_node_id, *node_id, atom.distance(&neighbor_atom));
        }
    }
    let save_fname = fname.replace(".pdb", "_graph.json");
    debug!("Parsing protein {}, node couunt {}. edge count {}", fname, protein_graph.node_count(), protein_graph.edge_count());

    let json = serde_json::to_string(&protein_graph)?;
    let mut file = File::create(&save_fname)?;
    debug!("Saved graph file {}", &save_fname);
    file.write_all(&json.as_bytes())?;

    Ok(())
}


fn main() {
    env_logger::init();
    let cmd = clap::Command::new("graphein")
        .bin_name("graphein")
        .arg(
            clap::arg!(--"pdb-glob" <PATH> "Glob pattern for protein files")
                .value_parser(clap::value_parser!(std::path::PathBuf)),
        )
        .arg(
            clap::arg!(--"cutoff" <f64> "Cutoff distance for graph edges")
                .value_parser(clap::value_parser!(f64)).default_value("3.5"),
        );
    

    let matches = cmd.get_matches();

    let edge_max_dist = matches.get_one::<f64>("cutoff").unwrap();
    let pdb_glob = glob(matches.get_one::<std::path::PathBuf>("pdb-glob").unwrap().to_str().unwrap()).expect("Failed to read glob pattern");

    let paths: Vec<String> = pdb_glob.map(|p| String::from(p.unwrap().to_str().unwrap())).collect();

    let results: Vec<Result<()>> = paths.par_iter().map(|p| process_pdb_file(p, edge_max_dist)).collect();

    let ok_res = results.iter().filter(|r| r.is_ok()).count();
    let err_res = results.iter().filter(|r| r.is_err()).count();

    info!("Processed {} proteins, failed {} times", ok_res, err_res);
    for e in results.iter().filter(|r| r.is_err()) {
        warn!("{:?}", e);
    }


    
}