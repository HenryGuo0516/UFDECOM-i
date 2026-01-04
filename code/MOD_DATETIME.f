	MODULE DATETIME1

	implicit none

c ====== New Types

	type Date
		integer Year
		integer Month
		integer Day
	endtype
	
	type Time
		integer Year
		integer Month
		integer Day
		integer Hour
		integer Minute
		integer Second
	endtype
	
c ====== Interfaces

	Interface operator(+)
		module procedure add_date_days
		module procedure add_days_date
		module procedure add_time_days0
		module procedure add_time_days1
		module procedure add_time_days2
		module procedure add_days_time0
		module procedure add_days_time1
		module procedure add_days_time2
		module procedure add_date_time
		module procedure add_time_date
		module procedure add_time_time
	endInterface

	Interface operator(-)
		module procedure minus_date_days
		module procedure minus_date_date
		module procedure minus_time_days0
		module procedure minus_time_days1
		module procedure minus_time_days2
		module procedure minus_time_time
		module procedure minus_time_date
		module procedure minus_date_time
	endInterface

c	Interface add_time_days
c		module procedure add_time_days0
c		module procedure add_time_days1
c		module procedure add_time_days2
c	endInterface

c	Interface add_days_time
c		module procedure add_days_time0
c		module procedure add_days_time1
c		module procedure add_days_time2
c	endInterface

c	Interface minus_time_days
c		module procedure minus_time_days0
c		module procedure minus_time_days1
c		module procedure minus_time_days2
c	endInterface

	Interface DaysToTime
		module procedure DaysToTime0
		module procedure DaysToTime1
		module procedure DaysToTime2
	endInterface

	Interface IsValid
		module procedure IsValidDate0
		module procedure IsValidTime0
		module procedure IsValidDate3I
		module procedure IsValidTime6I
	endInterface
	

	Interface IsValidDate
		module procedure IsValidDate0
		module procedure IsValidDate3I
	endInterface

	Interface IsValidTime
		module procedure IsValidTime0
		module procedure IsValidTime6I
	endInterface

	Interface NewTime
		module procedure NewTime6I
		module procedure NewTime3I
	endInterface

	Interface WeekDay
		module procedure WeekDay0
		module procedure WeekDay1
		module procedure WeekDay2
	endInterface
	
	contains

c ====== Functions
c --- Difine New Value:
	
	type(Date) function NewDate(Y,M,D)
		integer Y,M,D
		NewDate%Year=Y
		NewDate%Month=M
		NewDate%Day=D
		return
	endfunction

	type(Time) function NewTime6I(Y,Month,D,H,Minute,S)
		integer Y,Month,D,H,Minute,S
		NewTime6I%Year=Y
		NewTime6I%Month=Month
		NewTime6I%Day=D
		NewTime6I%Hour=H
		NewTime6I%Minute=Minute
		NewTime6I%Second=S
		return
	endfunction

	type(Time) function NewTime3I(H,Minute,S)
		integer H,Minute,S
		NewTime3I%Year=-1
		NewTime3I%Month=1
		NewTime3I%Day=1
		NewTime3I%Hour=H
		NewTime3I%Minute=Minute
		NewTime3I%Second=S
		return
	endfunction

c --- Logical Fuctions:

	Logical function IsLeapYear(Year)
		Integer Year,Y
		If (Year>0) Then
			Y=Year
		Else
			Y=Year+1
		Endif
		if (MODULO(Y,400)==0) then
			IsLeapYear=.True.
		elseif 	((MODULO(Y,100)/=0).AND.(MODULO(Y,4)==0)) then
			IsLeapYear=.True.
		else
			IsLeapYear=.False.
		endif
		return
	endfunction


c === IsValidDate

	Logical function IsValidDate0(Date1)
		type(Date) Date1
		integer Y,M,D
		Y=Date1%Year
		M=Date1%Month
		D=Date1%Day
		IsValidDate0=.True.
		If (M>12.OR.M<1.OR.D<1) IsValidDate0=.False.
		If (Y==0) IsValidDate0=.False.
		SelectCase (M)
		Case (1,3,5,7,8,10,12)
			If (D>31) IsValidDate0=.False.
		Case (4,6,9,11)
			If (D>30) IsValidDate0=.False.
		Case (2)
			If (IsLeapYear(Y)) then
				If (D>29) IsValidDate0=.False.
			Else
				If (D>28) IsValidDate0=.False.
			Endif
		EndSelect
		return
	endfunction

	Logical function IsValidDate3I(Y,M,D)
		Integer Y,M,D
		IsValidDate3I=IsValidDate0(NewDate(Y,M,D))
		return
	endfunction

c === IsValidTime

	Logical function IsValidTime0(Time1)
		type(time) Time1
		integer Y,Month,D,H,Minute,S
		Y=Time1%Year
		Month=Time1%Month
		D=Time1%Day
		H=Time1%Hour
		Minute=Time1%Minute
		S=Time1%Second
		IsValidTime0=.True.
		If (.Not.IsValidDate(Y,Month,D)) Then
			IsValidTime0=.False.
		Elseif (H<0.or.H>23.or.Minute<0.or.Minute>59.or.S<0.or.S>59)
     *         Then
			IsValidTime0=.False.
		Endif
		return
	endfunction

	Logical function IsValidTime6I(Y,Month,D,H,Minute,S)
		integer Y,Month,D,H,Minute,S
		IsValidTime6I=IsValidTime0(NewTime(Y,Month,D,H,Minute,S))
		return
	endfunction

c --- Convertions:

	type(time) function DateToTime(date1)
		type(date) date1
		DateToTime%Year=date1%Year
		DateToTime%Month=date1%Month
		DateToTime%Day=date1%Day
		DateToTime%Hour=0
		DateToTime%Minute=0
		DateToTime%Second=0
		return
	endfunction

	type(date) function TimeToDate(time1)
		type(time) time1
		TimeToDate%Year=time1%Year
		TimeToDate%Month=time1%Month
		TimeToDate%Day=time1%Day
		return
	endfunction

	Integer function DateToDays(Date1)
c从公元前1年开始的天数,采用简单定义的公历。
		type(Date) Date1
		Integer Y,Year,M,D,Days
		Integer a0,b0,a1,b1,a2,b2,a3,b3
		if (IsValid(Date1)) then
			Year=Date1%Year
			M=Date1%Month
			D=Date1%Day
			If (Year<0) Then
				Y=Year+1
			Else
				Y=Year
			Endif
			Select Case (M)
			Case (1)
				Days=0
			Case (2)
				Days=31
			Case (3)
				Days=59
			Case (4)
				Days=90
			Case (5)
				Days=120
			Case (6)
				Days=151
			Case (7)
				Days=181
			Case (8)
				Days=212
			Case (9)
				Days=243
			Case (10)
				Days=273
			Case (11)
				Days=304
			Case (12)
				Days=334
			End Select
			If (IsLeapYear(Year).and.M<=2) Days=Days-1 !减去闰年预加的1
			Days=Days+Y*365+floor(real(Y)/4.)-floor(real(Y)/100.)
     *+floor(real(Y)/400.)+D
c			天数=(月天数)+(年天数闰年预加1)+(日天数D-1)
			DateToDays=Days
		else
			DateToDays=0
		endif
		return
	endfunction DateToDays

	type(Date) function DaysToDate(Days) !从公元前1年开始,采用简单定义的公历。
		Integer Days,Y,Year,Month,Day
		Integer a0,b0,a1,b1,a2,b2,a3,b3
		Logical BC
		If (Days>=366) Then
			BC=.false.
		Else
			BC=.True.
		Endif
		b0=FLOOR(REAL(Days)/146097.)
		a0=modulo(Days,146097)
		if (a0<36525) then
			b1=0
			a1=a0
			b2=a1/1461
			a2=modulo(a1,1461)
			if (a2<366) then
				b3=0
				a3=a2
			else
				b3=(a2-1)/365
				a3=modulo(a2-1,365)
			endif
		else
			b1=(a0-1)/36524
			a1=modulo(a0-1,36524)
			if (a1<1460) then !前4年中
				b2=0
				a2=a1
				b3=a2/365
				a3=modulo(a2,365)
			else
				b2=(a1+1)/1461
				a2=modulo(a1+1,1461)
				if (a2<366) then
					b3=0
					a3=a2
				else
					b3=(a2-1)/365
					a3=modulo(a2-1,365)
				endif
			endif
		endif
		Y=b0*400+b1*100+b2*4+b3
		If (BC) Then
			Year=Y-1
		Else
			Year=Y
		Endif
c         --- Month Day:
		If (IsLeapYear(Year)) Then
			if (a3<31) then !1
				Month=1
				Day=a3+1
			elseif (a3<60) then !2
				Month=2
				Day=a3-30
			elseif (a3<91) then !3
				Month=3
				Day=a3-59
			elseif (a3<121) then !4
				Month=4
				Day=a3-90
			elseif (a3<152) then !5
				Month=5
				Day=a3-120
			elseif (a3<182) then !6
				Month=6
				Day=a3-151
			elseif (a3<213) then !7
				Month=7
				Day=a3-181
			elseif (a3<244) then !8
				Month=8
				Day=a3-212
			elseif (a3<274) then !9
				Month=9
				Day=a3-243
			elseif (a3<305) then !10
				Month=10
				Day=a3-273
			elseif (a3<335) then !11
				Month=11
				Day=a3-304
			else !12
				Month=12
				Day=a3-334
			endif
		Else
			if (a3<31) then !1
				Month=1
				Day=a3+1
			elseif (a3<59) then !2
				Month=2
				Day=a3-30
			elseif (a3<90) then !3
				Month=3
				Day=a3-58
			elseif (a3<120) then !4
				Month=4
				Day=a3-89
			elseif (a3<151) then !5
				Month=5
				Day=a3-119
			elseif (a3<181) then !6
				Month=6
				Day=a3-150
			elseif (a3<212) then !7
				Month=7
				Day=a3-180
			elseif (a3<243) then !8
				Month=8
				Day=a3-211
			elseif (a3<273) then !9
				Month=9
				Day=a3-242
			elseif (a3<304) then !10
				Month=10
				Day=a3-272
			elseif (a3<334) then !11
				Month=11
				Day=a3-303
			else !12
				Month=12
				Day=a3-333
			endif
		Endif
		DaysToDate%Year=Year
		DaysToDate%Month=Month
		DaysToDate%Day=Day
		return
	endfunction DaysToDate

	real(8) function TimeToDays(Time1)
		type(time) Time1
		integer Days
		integer Y,Month,D,H,Minute,S
		Y=Time1%Year
		Month=Time1%Month
		D=Time1%Day
		H=Time1%Hour
		Minute=Time1%Minute
		S=Time1%Second
		Days=DateToDays(NewDate(Y,Month,D))
		TimeToDays=DBLE(Days)
     *         +DBLE(H)/24._8+DBLE(Minute)/1440._8+DBLE(S)/86400._8 !release 2
		return
	endfunction

	type(time) function DaysToTime0(Days)
		real(8) Days
		type(date) Date1
		integer Y,Month,D,H,Minute,S,IDays
		real rH,rM,rS
		IDays=FLOOR(REAL(Days))
		Date1=DaysToDate(IDays)
		Y=Date1%Year
		Month=Date1%Month
		D=Date1%Day
		rH=(Days-IDays)*24.
		H=int(rH)
		rM=(rH-H)*60.
		Minute=int(rM)
		rS=(rM-Minute)*60.
		S=int(rS)
		DaysToTime0=NewTime(Y,Month,D,H,Minute,S)
		return
	endfunction

	type(time) function DaysToTime1(Days)
		real Days
		DaysToTime1=DaysToTime0(DBLE(Days))
		return
	endfunction

	type(time) function DaysToTime2(Days)
		integer Days
		DaysToTime2=DaysToTime0(DBLE(Days))
		return
	endfunction

c --- Operators:

c === date:

	type(date) function add_days_date(days1,date1)
		integer,intent(in) :: days1
		type(date),intent(in) :: date1
		add_days_date=DaysToDate(days1+DateToDays(date1))
		return
	endfunction

	type(date) function add_date_days(date1,days1)
		integer,intent(in) :: days1
		type(date),intent(in) :: date1
		add_date_days=DaysToDate(days1+DateToDays(date1))
		return
	endfunction

	type(date) function minus_date_days(date1,days1)
		integer,intent(in) :: days1
		type(date),intent(in) :: date1
		minus_date_days=DaysToDate(DateToDays(date1)-days1)
		return
	endfunction

	integer function minus_date_date(date1,date2)
		type(date),intent(in) :: date1,date2
		minus_date_date=DateToDays(date1)-DateToDays(date2)
		return
	endfunction

c === time:
	type(time) function add_date_time(date1,time2)
		type(time),intent(in) ::  time2
		type(date),intent(in) ::  date1
		add_date_time=
     *DaysToTime(TimeToDays(time2)+DBLE(DateToDays(date1)))	
		return
	endfunction

	type(time) function add_time_date(time2,date1)
		type(time),intent(in) ::  time2
		type(date),intent(in) ::  date1
		add_time_date=
     *DaysToTime(TimeToDays(time2)+DBLE(DateToDays(date1)))	
		return
	endfunction

	type(time) function add_time_time(time1,time2)
		type(time),intent(in) ::  time1,time2
		add_time_time=
     *DaysToTime(TimeToDays(time1)+TimeToDays(time2))	
		return
	endfunction

	type(time) function add_time_days0(time1,days1)
		type(time),intent(in) ::  time1
		real(8),intent(in) ::  days1
		add_time_days0=DaysToTime(TimeToDays(time1)+days1)
		return	
	endfunction

	type(time) function add_time_days1(time1,days1)
		type(time),intent(in) ::  time1
		real,intent(in) ::  days1	
		add_time_days1=DaysToTime(TimeToDays(time1)+DBLE(days1))
		return
	endfunction

	type(time) function add_time_days2(time1,days1)
		type(time),intent(in) ::  time1
		integer,intent(in) :: days1	
		add_time_days2=DaysToTime(TimeToDays(time1)+DBLE(days1))
		return
	endfunction

	type(time) function add_days_time0(days1,time1)
		type(time),intent(in) ::  time1
		real(8),intent(in) ::  days1
		add_days_time0=DaysToTime(TimeToDays(time1)+days1)
		return	
	endfunction

	type(time) function add_days_time1(days1,time1)
		type(time),intent(in) ::  time1
		real,intent(in) ::  days1	
		add_days_time1=DaysToTime(TimeToDays(time1)+DBLE(days1))
		return
	endfunction

	type(time) function add_days_time2(days1,time1)
		type(time),intent(in) ::  time1
		integer,intent(in) :: days1	
		add_days_time2=DaysToTime(TimeToDays(time1)+DBLE(days1))
		return
	endfunction

	type(time) function minus_time_days0(time1,days1)
		type(time),intent(in) ::  time1
		real(8),intent(in) ::  days1
		minus_time_days0=DaysToTime(TimeToDays(time1)-days1)
		return	
	endfunction

	type(time) function minus_time_days1(time1,days1)
		type(time),intent(in) ::  time1
		real,intent(in) ::  days1
		minus_time_days1=DaysToTime(TimeToDays(time1)-DBLE(days1))
		return	
	endfunction

	type(time) function minus_time_days2(time1,days1)
		type(time),intent(in) ::  time1
		integer,intent(in) ::  days1
		minus_time_days2=DaysToTime(TimeToDays(time1)-DBLE(days1))
		return	
	endfunction

	real(8) function minus_time_time(time1,time2)
		type(time),intent(in) ::  time1,time2
		minus_time_time=TimeToDays(time1)-TimeToDays(time2)
		return
	endfunction	

	real(8) function minus_time_date(time1,date1)
		type(time),intent(in) ::  time1
		type(date),intent(in) ::  date1
		minus_time_date=TimeToDays(time1)-DBLE(DateToDays(date1))
		return
	endfunction	

	real(8) function minus_date_time(date1,time1)
		type(time),intent(in) ::  time1
		type(date),intent(in) ::  date1
		minus_date_time=DBLE(DateToDays(date1))-TimeToDays(time1)
		return
	endfunction	

c === WeakDay:
	integer function WeekDay0(Date1)
		type(date),intent(in) :: date1
		WeekDay0=modulo(DateToDays(date1)-1,7)
		return
	endfunction	

	integer function WeekDay1(Time1)
		type(time),intent(in) ::  Time1
		WeekDay1=WeekDay0(TimeToDate(Time1))
		return
	endfunction	

	integer function WeekDay2(Y,M,D)
		integer,intent(in) :: Y,M,D
		If (IsValidDate(Y,M,D)) Then
			WeekDay2=WeekDay0(NewDate(Y,M,D))
		Else
			WeekDay2=-1
		Endif	
		return
	endfunction	


	ENDMODULE DATETIME1
