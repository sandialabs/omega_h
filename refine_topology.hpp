void refine_domains_to_pairs(
    Mesh& mesh,
    Int dim,
    LOs keys2edges,
    LOs keys2midverts,
    LOs old_verts2new_verts,
    LOs& keys2pairs,
    LOs& pair_verts2verts);

void refine_domains_to_cuts(
    Mesh& mesh,
    Int dim,
    LOs keys2edges,
    LOs keys2midverts,
    LOs old_verts2new_verts,
    LOs& keys2cuts,
    LOs& cut_verts2verts);

void combine_pairs_and_cuts(
    Int ent_dim,
    LOs keys2cuts,
    LOs keys2pairs,
    LOs cut_verts2verts,
    LOs pair_verts2verts,
    LOs& keys2prods,
    LOs& prod_verts2verts);

void refine_products(
    Mesh& mesh,
    Int ent_dim,
    LOs keys2edges,
    LOs keys2midverts,
    LOs old_verts2new_verts,
    LOs& keys2prods,
    LOs& prod_verts2verts);
