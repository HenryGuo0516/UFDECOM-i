#include "DEFS.h"
	SUBROUTINE CAL_TH

	USE MOD_GLOBAL
      CHARACTER*10 YEAR_S,MONTH_S,DAY_S,HOUR_S,MINUTE_S,SECOND_S
      DELTA_TGTH_S=DBLE((TG_TH-TDAY-2)*24.*60.*60.)
      DELTA_TGTH_S1=DBLE((TG_TH-TDAY)*24.*60.*60.)
      DELTA_TGTH_S7=DBLE((TG_TH-TDAY+5)*24.*60.*60.)
      IF (OUT_T0.EQ..FALSE..AND.ABS(DELTA_TGTH_S1)<=5.*DTI) THEN
          OPEN(10024,FILE='TONGHUA.TXT')
          WRITE(10024,*) NSTEP,NN_FPT,TDAY
          CLOSE(10024)
          CALL NMC_DIFF
          CALL DENS
          WRITE(YEAR_S  ,'(I4.4)')YEAR_TH
          WRITE(MONTH_S ,'(I2.2)')MONTH_TH
          WRITE(DAY_S   ,'(I2.2)')DAY_TH
          WRITE(HOUR_S  ,'(I2.2)')HOUR_TH
          WRITE(MINUTE_S,'(I2.2)')MINUTE_TH
          WRITE(SECOND_S,'(I2.2)')SECOND_TH
          FN='restart_'//TRIM(YEAR_S)//TRIM(MONTH_S)
     *        //TRIM(DAY_S)//TRIM(HOUR_S)//TRIM(MINUTE_S)
     *        //TRIM(SECOND_S)
#include "WRITE_RST.h"
          OUT_T0=.TRUE.
      ELSEIF (OUT_T07.EQ..FALSE..AND.ABS(DELTA_TGTH_S7)<=5.*DTI) THEN
          WRITE(YEAR_S  ,'(I4.4)')YEAR_TH
          WRITE(MONTH_S ,'(I2.2)')MONTH_TH
          WRITE(DAY_S   ,'(I2.2)')DAY_TH
          WRITE(HOUR_S  ,'(I2.2)')HOUR_TH
          FN='restart_P1_'//TRIM(YEAR_S)//TRIM(MONTH_S)
     *        //TRIM(DAY_S)//TRIM(HOUR_S)
#include "WRITE_RST.h"
          OUT_T07=.TRUE.
      ELSEIF (OUT_T01.EQ..FALSE..AND.ABS(DELTA_TGTH_S)<=5.*DTI) THEN
          WRITE(YEAR_S  ,'(I4.4)')YEAR_TH
          WRITE(MONTH_S ,'(I2.2)')MONTH_TH
          WRITE(DAY_S   ,'(I2.2)')DAY_TH
          WRITE(HOUR_S  ,'(I2.2)')HOUR_TH
          FN='restart_M2_'//TRIM(YEAR_S)//TRIM(MONTH_S)
     *        //TRIM(DAY_S)//TRIM(HOUR_S)
#include "WRITE_RST.h"
          OUT_T01=.TRUE.
      ENDIF
	RETURN  
	END SUBROUTINE CAL_TH