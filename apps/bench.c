#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libs/cLOG/cLOG.h"

#include "core/CFD.h"

char fileNames[][50] = {
    "40_40_100_UDS_SECOND_07_07",
    "40_40_400_UDS_SECOND_07_07",
    "40_40_1000_UDS_SECOND_07_07",
    "40_40_1000_CDS_SECOND_07_07",
    "40_40_1000_QUICK_SECOND_07_07",
    "80_80_1000_UDS_SECOND_07_07",
    "80_80_1000_CDS_SECOND_07_07",
    "80_80_1000_QUICK_SECOND_07_07",
    "129_129_1000_UDS_SECOND_07_07",
    "129_129_1000_CDS_SECOND_07_07",
    "129_129_1000_QUICK_SECOND_07_07",
    "129_129_1000_UDS_FOURTH_07_07",
    "129_129_1000_CDS_FOURTH_07_07",
    "129_129_1000_QUICK_FOURTH_07_07",
};

int skipping[] = {
    -1,
    // 0,
    // 1,
    // 2,
    // 3,
    // 4,
    // 5,
    // 6,
    // 7,
    // 8,
    // 9,
    // 10,
    // 11,
    // 12,
    // 13,
};

bool isSkipping(uint8_t index);

int main(int argc, char *argv[])
{
    log_set_level(LOG_INFO);

    for (uint8_t i = 0; i < sizeof(fileNames) / sizeof(fileNames[0]); i++)
    {

        if (isSkipping(i))
        {
            log_info("#############################################");
            log_info("Simulation skipping %02d_%s", i, fileNames[i]);
            log_info("#############################################");
            log_info("\n\n");

            continue;
        }

        CFD_t *cfd = CFD_Allocate();

        log_info("#############################################");
        log_info("Simulation starting (%02d)", i);
        log_info("#############################################");

        sprintf(cfd->in->file->path, "%s", DEFAULT_IN_FILE_PATH);
        sprintf(cfd->in->file->name, "%02d_%s", i, fileNames[i]);
        cfd->in->file->extension = FILE_String_to_Extension(DEFAULT_IN_FILE_FORMAT);

        CFD_Prepare(cfd, argc, argv);

        CFD_Solve(cfd);

        CFD_Finalize(cfd);

        log_info("#############################################");
        log_info("Simulation ended (%02d)", i);
        log_info("#############################################");
        log_info("\n\n");

        CFD_Free(cfd);
    }

    return EXIT_SUCCESS;
}

bool isSkipping(uint8_t index)
{
    for (size_t i = 0; i < sizeof(skipping) / sizeof(skipping[0]); i++)
    {
        if (index == skipping[i])
        {
            return true;
        }
    }
    return false;
}
