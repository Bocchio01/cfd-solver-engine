typedef struct CFD_t CFD_t;

#include "../SIMPLE.h"
#include "U_equation.h"
#include "../../../schemes/schemes.h"

#include "libs/cALGEBRA/cMAT2D.h"
#include "../../../../CFD.h"

#include <stdint.h>

void CFD_SIMPLE_U_Compute_Coefficients(CFD_t *cfd, SIMPLE_t *simple)
{
    uint16_t i;
    uint16_t j;

    for (cfd->engine->method->index->j = cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j++)
    {
        for (cfd->engine->method->index->i = cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1;
             cfd->engine->method->index->i++)
        {
            i = cfd->engine->method->index->i;
            j = cfd->engine->method->index->j;

            // Diagonal of the implicit operator: UDS convection + 2nd order diffusion. It is always dominant, while the
            // diagonal of QUICK/CDS can be small or even negative at high cell Peclet numbers: QUICK/CDS enter only
            // through the neighbour coefficients of the residual (deferred correction).
            CFD_Scheme_Convection_UDS(cfd, i, j, u);
            CFD_Scheme_Diffusion_SECOND(cfd);
            simple->diagonal->u->data[j][i] = cfd->engine->schemes->convection->coefficients->data[PP] + cfd->engine->schemes->diffusion->coefficients->data[PP];

            if (i == cfd->engine->mesh->n_ghosts ||
                j == cfd->engine->mesh->n_ghosts ||
                i == cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1 ||
                j == cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts - 1)
            {
                CFD_Scheme_Convection_UDS(cfd, i, j, u);
                CFD_Scheme_Diffusion_SECOND(cfd);
            }
            else
            {
                cfd->engine->schemes->convection->callable(cfd, i, j, u);
                cfd->engine->schemes->diffusion->callable(cfd);
            }

            for (uint16_t k = 0; k < cfd->engine->schemes->convection->coefficients->length; k++)
            {
                simple->Ap_coefficients->u->data[j][i][k] = cfd->engine->schemes->convection->coefficients->data[k] + cfd->engine->schemes->diffusion->coefficients->data[k];
            }
        }
    }
}

void CFD_SIMPLE_U_Apply_BC(CFD_t *cfd)
{
    // Horizontal walls
    for (uint16_t i = cfd->engine->mesh->n_ghosts;
         i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts;
         i++)
    {
        cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts - 1][i] = 2.0 * 0.0 - cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts][i];
        cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Ny][i] = 2.0 * cfd->in->uLid - cfd->engine->method->state->u->data[cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Ny - 1][i];
    }

    // Vertical walls
    for (uint16_t j = cfd->engine->mesh->n_ghosts;
         j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
         j++)
    {
        cfd->engine->method->state->u->data[j][cfd->engine->mesh->n_ghosts - 1] = 0.0;
        cfd->engine->method->state->u->data[j][cfd->engine->mesh->n_ghosts + cfd->engine->mesh->nodes->Nx - 1] = 0.0;
    }
}

void CFD_SIMPLE_U_Compute_Residuals(CFD_t *cfd, SIMPLE_t *simple)
{
    uint16_t i;
    uint16_t j;
    double phi;

    for (cfd->engine->method->index->j = cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j++)
    {
        for (cfd->engine->method->index->i = cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1;
             cfd->engine->method->index->i++)
        {
            i = cfd->engine->method->index->i;
            j = cfd->engine->method->index->j;
            simple->residual->u->data[j][i] = 0.0;

            for (uint16_t k = 0; k < simple->Ap_coefficients->u->depth; k++)
            {
                phi = CFD_Get_State(cfd, u, i + ((k / 5) % 5) - 2, j + (k % 5) - 2);
                if (k == PP)
                {
                    simple->residual->u->data[j][i] -= simple->Ap_coefficients->u->data[j][i][k] * phi;
                }
                else
                {
                    simple->residual->u->data[j][i] += simple->Ap_coefficients->u->data[j][i][k] * phi;
                }
            }

            simple->residual->u->data[j][i] += (CFD_Get_State(cfd, p, i, j) - CFD_Get_State(cfd, p, i + 1, j)) * cfd->engine->mesh->element->size->dy;
        }
    }
}

void CFD_SIMPLE_U_Compute_State(CFD_t *cfd, SIMPLE_t *simple)
{
    uint16_t i;
    uint16_t j;
    double sum;

    // u of the previous outer iteration: it is the reference of the under-relaxation for every sweep.
    // Relaxing towards the value of the previous sweep would make the relaxation fade as the sweeps accumulate.
    MAT2D_Copy(cfd->engine->method->state->u, simple->state->u);

    for (uint8_t sweep = 0; sweep < simple->number_of_sweeps->u; sweep++)
    {

        for (cfd->engine->method->index->j = cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->j++)
        {
            for (cfd->engine->method->index->i = cfd->engine->mesh->n_ghosts;
                 cfd->engine->method->index->i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1;
                 cfd->engine->method->index->i++)
            {
                i = cfd->engine->method->index->i;
                j = cfd->engine->method->index->j;

                sum = 0.0;

                for (uint8_t k = 0; k < simple->Ap_coefficients->u->depth; k++)
                {
                    if (k != PP)
                    {
                        sum += simple->Ap_coefficients->u->data[j][i][k] * CFD_Get_State(cfd, u, i + ((k / 5) % 5) - 2, j + (k % 5) - 2);
                    }
                }

                sum += (CFD_Get_State(cfd, p, i, j) - CFD_Get_State(cfd, p, i + 1, j)) * cfd->engine->mesh->element->size->dy;

                // Gauss-Seidel update with the diagonal of the implicit operator (for UDS it is A_P, so this is sum / A_P)
                sum -= simple->Ap_coefficients->u->data[j][i][PP] * CFD_Get_State(cfd, u, i, j);
                cfd->engine->method->state->u->data[j][i] = simple->under_relaxation->u * (CFD_Get_State(cfd, u, i, j) + sum / simple->diagonal->u->data[j][i]) + (1.0 - simple->under_relaxation->u) * simple->state->u->data[j][i];
            }
        }
    }
}

void CFD_SIMPLE_U_Compute_Correct_States(CFD_t *cfd, SIMPLE_t *simple)
{
    uint16_t i;
    uint16_t j;
    double u_prime;

    for (cfd->engine->method->index->j = cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
         cfd->engine->method->index->j++)
    {
        for (cfd->engine->method->index->i = cfd->engine->mesh->n_ghosts;
             cfd->engine->method->index->i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts - 1;
             cfd->engine->method->index->i++)
        {
            i = cfd->engine->method->index->i;
            j = cfd->engine->method->index->j;

            u_prime = cfd->engine->mesh->element->size->dy / simple->diagonal->u->data[j][i] * (simple->state->pp->data[j][i] - simple->state->pp->data[j][i + 1]);
            cfd->engine->method->state->u->data[j][i] = CFD_Get_State(cfd, u, i, j) + simple->under_relaxation->u * u_prime;
        }
    }
}