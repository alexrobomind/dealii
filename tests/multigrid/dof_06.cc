#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <vector>

constexpr int DIM = 2;

using namespace dealii;

template <typename FaceType>
void
print_face(const FaceType &face, unsigned int level)
{
  AssertThrow(
    face->get_fe(0).n_dofs_per_face() > 0,
    ExcMessage(
      "Please use a finite element with at least 1 DoF for this test"));

  std::vector<types::global_dof_index> mg_dofs(
    face->get_fe(0).n_dofs_per_face());
  face->get_mg_dof_indices(level, mg_dofs);

  // Print face DoFs
  std::cout << mg_dofs[0];
  for (unsigned int i_dof = 1; i_dof < mg_dofs.size(); ++i_dof)
    std::cout << ", " << mg_dofs[i_dof];
}

template <int dim>
void
check()
{
  typedef Triangulation<dim> Tria;
  // Create triangulation
  Tria tria(Tria::limit_level_difference_at_vertices);

  // Create coarse grid
  GridGenerator::hyper_cube(tria, -1, 1, true);

  // Create periodic boundary along x axis
  std::vector<GridTools::PeriodicFacePair<typename Tria::cell_iterator>>
    face_pairs;

  GridTools::collect_periodic_faces(tria, 0, 1, 0, face_pairs);

  tria.add_periodicity(face_pairs);

  // Refine once and then refine the corner again
  tria.refine_global(1);

  tria.begin_active()->set_refine_flag();
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();

  // Assemble DoFs
  DoFHandler<dim> dof_handler(tria);
  FE_Q<dim>       fe_q(1);

  dof_handler.distribute_dofs(fe_q);
  dof_handler.distribute_mg_dofs();

  // Try to assemble MGConstrainedDofs
  MGConstrainedDoFs constrained_dofs;
  constrained_dofs.initialize(dof_handler);

  // Report which DIM is ok
  std::cout << "--- DIM " << dim << " OK ---" << std::endl;

  // Debug printouts:
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      std::cout << "Level " << level << ": " << std::endl;

      // Print faces that connect to coarser DoF
      std::cout << "  Faces neighboring coarser cells:" << std::endl;
      for (typename DoFHandler<dim>::cell_iterator level_cell =
             dof_handler.begin(level);
           level_cell != dof_handler.end(level);
           ++level_cell)
        {
          // Iterate over all faces
          for (unsigned int i_face = 0; i_face < 2 * dim; ++i_face)
            {
              if ((level_cell->at_boundary(i_face) &&
                   !level_cell->has_periodic_neighbor(i_face)) ||
                  level_cell->neighbor_or_periodic_neighbor(i_face)->level() ==
                    level_cell->level())
                continue;

              // Get face
              auto face = level_cell->face(i_face);

              std::cout << "    ";
              print_face(face, level);

              // Print whether the neighbor is periodic
              if (level_cell->has_periodic_neighbor(i_face))
                std::cout << " via periodic boundary" << std::endl;
              else
                std::cout << " via internal face" << std::endl;
            }
        }

      // Print hanging node DoFs
      std::cout << "  Refinement edge DoFs: ";
      constrained_dofs.get_refinement_edge_indices(level).print(std::cout);
      std::cout << std::endl;

      // Print faces that have periodic neighbors
      std::cout << "  Faces having periodic neighbors of same level:"
                << std::endl;
      for (typename DoFHandler<dim>::cell_iterator level_cell =
             dof_handler.begin(level);
           level_cell != dof_handler.end(level);
           ++level_cell)
        {
          // Iterate over all faces
          for (unsigned int i_face = 0; i_face < 2 * dim; ++i_face)
            {
              if (!level_cell->has_periodic_neighbor(i_face) ||
                  level_cell->periodic_neighbor(i_face)->level() !=
                    level_cell->level())
                continue;

              typename DoFHandler<dim>::cell_iterator neighbor_cell =
                level_cell->periodic_neighbor(i_face);

              // Each side only once
              if (level_cell->index() > neighbor_cell->index())
                continue;

              unsigned int neighbor_i_face =
                level_cell->periodic_neighbor_face_no(i_face);

              // Get face
              auto face          = level_cell->face(i_face);
              auto neighbor_face = neighbor_cell->face(neighbor_i_face);

              std::cout << "    ";
              print_face(face, level);
              std::cout << " <-> ";
              print_face(neighbor_face, level);

              // Print whether the neighbor is periodic
              if (level_cell->has_periodic_neighbor(i_face))
                std::cout << " via periodic boundary" << std::endl;
              else
                std::cout << " via internal face" << std::endl;
            }
        }

      // Print constraint matrix
      std::cout << "  Constraint matrix:" << std::endl;
      constrained_dofs.get_level_constraint_matrix(level).print(std::cout);
    }
}

int
main()
{
  // check<1>();
  check<2>();
  check<3>();
  return 0;
}
