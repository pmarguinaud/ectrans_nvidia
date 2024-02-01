! (C) Copyright 2015- ECMWF.
! (C) Copyright 2015- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE VORDIV_TO_UV(PSPVOR,PSPDIV,PSPU,PSPV,KSMAX,KVSETUV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

REAL(KIND=JPRB), INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB), INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPU(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPV(:,:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)

CALL ABORT_TRANS('VORDIV_TO_UV IS NOT IMPLEMENTED!')

END SUBROUTINE VORDIV_TO_UV

