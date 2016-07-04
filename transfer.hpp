void transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh, LOs key_verts2verts,
                            LOs keys2midverts, LOs same_verts2old_verts,
                            LOs same_verts2new_verts);

void transfer_metric(Mesh* old_mesh, Mesh* new_mesh, LOs key_verts2verts,
                     LOs keys2midverts, LOs same_verts2old_verts,
                     LOs same_verts2new_verts);

void transfer_inherit_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
                             Int prod_dim, LOs keys2prods, LOs prods2new_ents,
                             LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_length(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
                     LOs same_ents2new_ents, LOs prods2new_ents);

void transfer_quality(Mesh* old_mesh, Mesh* new_mesh, LOs same_ents2old_ents,
                      LOs same_ents2new_ents, LOs prods2new_ents);

void transfer_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
                     LOs keys2midverts, Int prod_dim, LOs keys2prods,
                     LOs prods2new_ents, LOs same_ents2old_ents,
                     LOs same_ents2new_ents);

void transfer_no_products(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
                          LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_conserve(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
                       LOs keys2kds, LOs keys2prods, LOs prods2new_ents,
                       LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_coarsen(Mesh* old_mesh, Mesh* new_mesh, LOs keys2verts,
                      Adj keys2doms, Int prod_dim, LOs prods2new_ents,
                      LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_copy(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim);

void transfer_inherit_swap(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
                           LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
                           LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_swap(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim, LOs keys2edges,
                   LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
                   LOs same_ents2new_ents);
