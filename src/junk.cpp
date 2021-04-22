std::string Mesh::string() {
  std::ostringstream oss;
  oss << "Mesh:" 
      << "\ncomm()->size              = " << comm()->size()
      << "\nparting                   = " << parting()
      << "\ndim                       = " << dim()
      << "\nfamily                    = " << family()
      << "\nnents                     = " << nents(dim())
      << "\nnents(0)                  = " << nents(0)
      << "\nnelems                    = " << nelems()
      << "\nnregions                  = " << nregions()
      << "\nnfaces                    = " << nfaces()
      << "\nnedges                    = " << nedges()
      << "\nnverts                    = " << nverts()
      << "\nnglobal_ents              = " << nglobal_ents(dim())
      << "\nnglobal_ents(0)           = " << nglobal_ents(0)
      << "\nntags                     = " << ntags(dim())
      << "\nntags(0)                  = " << ntags(0)
      << "\nmin_quality               = " << min_quality()
      << "\nmax_length                = " << max_length()
      << "\ncould_be_shared           = " << could_be_shared(dim())
      << "\ncould_be_shared(0)        = " << could_be_shared(0)
      << "\nowners_have_all_upward    = " << owners_have_all_upward(dim())
      << "\nowners_have_all_upward(0) = " << owners_have_all_upward(0)
      << "\nhave_all_upward           = " << have_all_upward()
      << "\nimbalance                 = " << imbalance();
  return oss.str();
}

