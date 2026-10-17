#include "SIMPLE.h"

typedef struct CFD_t CFD_t;

#include <stdlib.h>
#include <math.h>
#include "../methods.h"
#include "../../../CFD.h"
#include "../../schemes/schemes.h"
#include "libs/cALGEBRA/cMAT3D.h"
#include "libs/cALGEBRA/cMAT2D.h"
#include "libs/cALGEBRA/cVEC.h"
#include "libs/cLOG/cLOG.h"
#include "P_equation.h"
#include "V_equation.h"
#include "U_equation.h"

void CFD_SIMPLE(CFD_t *cfd, cJSON *args)
{
    SIMPLE_t *simple;

    simple = CFD_SIMPLE_Allocate(cfd, args);

    for (cfd->engine->method->iteractions = 0;
         cfd->engine->method->iteractions < cfd->engine->method->maxIter;
         cfd->engine->method->iteractions++)
    {

        CFD_SIMPLE_Compute_U(cfd, simple);
        CFD_SIMPLE_Compute_V(cfd, simple);
        CFD_SIMPLE_Compute_P(cfd, simple);

        CFD_SIMPLE_U_Compute_Correct_States(cfd, simple);
        CFD_SIMPLE_V_Compute_Correct_States(cfd, simple);
        CFD_SIMPLE_P_Compute_Correct_States(cfd, simple);

        // MAT2D_Print_States(cfd->engine->method->state->u);
        // MAT2D_Print_States(cfd->engine->method->state->v);
        // MAT2D_Print_States(cfd->engine->method->state->p);

        // Convergence monitor (same norm as SCGS and as Eq. 161 of the report):
        // sum over the whole domain of the absolute value of the (signed) residuals.
        double residual_u = 0.0;
        double residual_v = 0.0;
        double residual_p = 0.0;

        for (uint16_t j = cfd->engine->mesh->n_ghosts;
             j < cfd->engine->mesh->nodes->Ny + cfd->engine->mesh->n_ghosts;
             j++)
        {
            for (uint16_t i = cfd->engine->mesh->n_ghosts;
                 i < cfd->engine->mesh->nodes->Nx + cfd->engine->mesh->n_ghosts;
                 i++)
            {
                residual_u += fabs(simple->residual->u->data[j][i]);
                residual_v += fabs(simple->residual->v->data[j][i]);
                residual_p += fabs(simple->residual->p->data[j][i]);
            }
        }

        // fmax() silently discards a NaN, so a diverged field would look "converged" (residual = 0):
        // the NaN is propagated explicitly.
        if (isnan(residual_u) || isnan(residual_v) || isnan(residual_p))
        {
            cfd->engine->method->residual->data[cfd->engine->method->iteractions] = NAN;
        }
        else
        {
            cfd->engine->method->residual->data[cfd->engine->method->iteractions] = fmax(fmax(residual_u, residual_v), residual_p);
        }

        // printf("Residual: %.10f\n", cfd->engine->method->residual->data[cfd->engine->method->iteractions]);

        if (isnan(cfd->engine->method->residual->data[cfd->engine->method->iteractions]) ||
            isinf(cfd->engine->method->residual->data[cfd->engine->method->iteractions]))
        {
            log_fatal("Algorithm diverged at iteration %d", cfd->engine->method->iteractions);
            CFD_SIMPLE_Free(simple);
            exit(EXIT_FAILURE);
        }

        if (cfd->engine->method->iteractions % 100 == 0 ||
            cfd->engine->method->residual->data[cfd->engine->method->iteractions] < cfd->engine->method->tolerance)
        {
            log_info("\nIteration:\t%d\nResidual:\t%.10f",
                     cfd->engine->method->iteractions,
                     cfd->engine->method->residual->data[cfd->engine->method->iteractions]);
        }

        if (cfd->engine->method->residual->data[cfd->engine->method->iteractions] < cfd->engine->method->tolerance &&
            cfd->engine->method->iteractions > 1)
        {
            log_info("Algorithm converged in %d iterations", cfd->engine->method->iteractions);
            CFD_SIMPLE_Free(simple);
            break;
        }
    }
}

SIMPLE_t *CFD_SIMPLE_Allocate(CFD_t *cfd, cJSON *args)
{
    SIMPLE_t *simple;
    uint16_t rows = cfd->engine->mesh->data->x->rows;
    uint16_t cols = cfd->engine->mesh->data->x->cols;

    simple = (SIMPLE_t *)malloc(sizeof(SIMPLE_t));

    if (simple != NULL)
    {
        simple->residual = (SIMPLE_residual_t *)malloc(sizeof(SIMPLE_residual_t));
        simple->state = (SIMPLE_state_t *)malloc(sizeof(SIMPLE_state_t));
        simple->diagonal = (SIMPLE_diagonal_t *)malloc(sizeof(SIMPLE_diagonal_t));
        simple->Ap_coefficients = (SIMPLE_Ap_coefficients_t *)malloc(sizeof(SIMPLE_Ap_coefficients_t));
        simple->number_of_sweeps = (SIMPLE_number_of_sweeps_t *)malloc(sizeof(SIMPLE_number_of_sweeps_t));
        simple->under_relaxation = (SIMPLE_under_relaxation_t *)malloc(sizeof(SIMPLE_under_relaxation_t));

        if (simple->residual != NULL &&
            simple->state != NULL &&
            simple->diagonal != NULL &&
            simple->Ap_coefficients != NULL &&
            simple->number_of_sweeps != NULL &&
            simple->under_relaxation != NULL)
        {
            simple->state->u = MAT2D_Init(rows, cols);
            simple->state->v = MAT2D_Init(rows, cols);
            simple->state->pp = MAT2D_Init(rows, cols);

            simple->diagonal->u = MAT2D_Init(rows, cols);
            simple->diagonal->v = MAT2D_Init(rows, cols);

            simple->residual->u = MAT2D_Init(rows, cols);
            simple->residual->v = MAT2D_Init(rows, cols);
            simple->residual->p = MAT2D_Init(rows, cols);

            simple->Ap_coefficients->u = MAT3D_Init(rows, cols, EENN + 1);
            simple->Ap_coefficients->v = MAT3D_Init(rows, cols, EENN + 1);
            simple->Ap_coefficients->pp = MAT3D_Init(rows, cols, EENN + 1);

            const cJSON *under_relaxation = cJSON_GetObjectItemCaseSensitive(args, "under_relaxation");
            const cJSON *u_urf = cJSON_GetObjectItemCaseSensitive(under_relaxation, "u");
            const cJSON *v_urf = cJSON_GetObjectItemCaseSensitive(under_relaxation, "v");
            const cJSON *p_urf = cJSON_GetObjectItemCaseSensitive(under_relaxation, "p");
            simple->under_relaxation->u = (float)(cJSON_IsNumber(u_urf) ? cJSON_GetNumberValue(u_urf) : DEFAULT_ENGINE_METHOD_UNDER_RELAXATION_U);
            simple->under_relaxation->v = (float)(cJSON_IsNumber(v_urf) ? cJSON_GetNumberValue(v_urf) : DEFAULT_ENGINE_METHOD_UNDER_RELAXATION_V);
            simple->under_relaxation->p = (float)(cJSON_IsNumber(p_urf) ? cJSON_GetNumberValue(p_urf) : DEFAULT_ENGINE_METHOD_UNDER_RELAXATION_P);

            const cJSON *sweeps = cJSON_GetObjectItemCaseSensitive(args, "number_of_sweeps");
            if (!cJSON_IsObject(sweeps))
            {
                sweeps = cJSON_GetObjectItemCaseSensitive(args, "sweeps");
            }
            const cJSON *u_sweep = cJSON_GetObjectItemCaseSensitive(sweeps, "u");
            const cJSON *v_sweep = cJSON_GetObjectItemCaseSensitive(sweeps, "v");
            const cJSON *pp_sweep = cJSON_GetObjectItemCaseSensitive(sweeps, "p");
            simple->number_of_sweeps->u = (float)(cJSON_IsNumber(u_sweep) ? cJSON_GetNumberValue(u_sweep) : DEFAULT_ENGINE_METHOD_SWEEPS_U);
            simple->number_of_sweeps->v = (float)(cJSON_IsNumber(v_sweep) ? cJSON_GetNumberValue(v_sweep) : DEFAULT_ENGINE_METHOD_SWEEPS_V);
            simple->number_of_sweeps->pp = (float)(cJSON_IsNumber(pp_sweep) ? cJSON_GetNumberValue(pp_sweep) : DEFAULT_ENGINE_METHOD_SWEEPS_P);

            return simple;
        }
    }

    log_fatal("Failed to allocate memory for SIMPLE method");
    exit(EXIT_FAILURE);
}

void CFD_SIMPLE_Free(SIMPLE_t *simple)
{
    MAT2D_Free(simple->state->u);
    MAT2D_Free(simple->state->v);
    MAT2D_Free(simple->state->pp);
    MAT2D_Free(simple->diagonal->u);
    MAT2D_Free(simple->diagonal->v);
    free(simple->residual);
    free(simple->state);
    free(simple->diagonal);
    free(simple);
}

void CFD_SIMPLE_Compute_U(CFD_t *cfd, SIMPLE_t *simple)
{
    CFD_SIMPLE_U_Compute_Coefficients(cfd, simple);
    CFD_SIMPLE_U_Apply_BC(cfd);
    CFD_SIMPLE_U_Compute_Residuals(cfd, simple);
    CFD_SIMPLE_U_Compute_State(cfd, simple);
}

void CFD_SIMPLE_Compute_V(CFD_t *cfd, SIMPLE_t *simple)
{
    CFD_SIMPLE_V_Compute_Coefficients(cfd, simple);
    CFD_SIMPLE_V_Apply_BC(cfd);
    CFD_SIMPLE_V_Compute_Residuals(cfd, simple);
    CFD_SIMPLE_V_Compute_State(cfd, simple);
}

void CFD_SIMPLE_Compute_P(CFD_t *cfd, SIMPLE_t *simple)
{
    CFD_SIMPLE_P_Compute_Coefficients(cfd, simple);
    CFD_SIMPLE_P_Apply_BC(cfd, simple);
    CFD_SIMPLE_P_Compute_Residuals(cfd, simple);
    CFD_SIMPLE_P_Compute_State(cfd, simple);
}