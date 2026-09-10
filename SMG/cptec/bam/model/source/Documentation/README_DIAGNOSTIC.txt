=========================================================================================

Instrucoes para anexar, outras variaveis diagnosticas no Modelo Global 
saulo@ime.usp.br
=========================================================================================

------------------------------------------------------
1 A Estrutura de diretorio do modelo Global deve ser:
------------------------------------------------------


                                         Global
                                           |
                                           |
                  ------------------------------------------------------
                  |                        |                           |
                model                     run                          pos
                  |                                                    |
       -----------------------                              ------------------------
       |       |      |       |                             |        |      |      |
    datain  dataout  exec  source                         datain  dataout  exec  source
               |                                                     |  
       ----------------                                       ---------------- 
       |              |                                       |              |    
    T042L28        T341L64                                 T042L28        T341L64
    
=========================================================================================
------------------------------------------------------
2 Module Constants:
------------------------------------------------------
  2.1 No Module Constants deve-se alterar o valor parameter ndavl
      este parameter indica quantas variaveis disgnosticas estao
      disponiveis.
      
      Por exemplo:
       ----------------- Modulo Constants -----------------
      
        INTEGER (KIND=i8), PARAMETER :: ndavl = 95  ! Number of Available Diagnostics
                                                    ! Fields

       ----------------- Modulo Constants -----------------
=========================================================================================
------------------------------------------------------
3 Module Diagnostics:
------------------------------------------------------

  3.1 O modulo Diagnostics se encontra no diretorio:
      " ~/Global/model/source/Diagnostics.f90"
      
      E necessario acrescentar um novo "parameter para" cada variavel diagnostica
      a ser inserida na tabela dos "Available Diagnostics Indexes"
      Por exemplo:
       ----------------- Modulo Diagnostics -----------------
        INTEGER, PUBLIC, PARAMETER :: nDiag_omegav =  7 ! omega
        INTEGER, PUBLIC, PARAMETER :: nDiag_sigdot =  8 ! sigma dot
        INTEGER, PUBLIC, PARAMETER :: nDiag_toprec =  9 ! total precipiation
        INTEGER, PUBLIC, PARAMETER :: nDiag_cvprec = 10 ! convective precipitation
        INTEGER, PUBLIC, PARAMETER :: nDiag_lsprec = 11 ! large scale precipitation 
        INTEGER, PUBLIC, PARAMETER :: nDiag_snowfl = 12 ! snowfall
        INTEGER, PUBLIC, PARAMETER :: nDiag_runoff = 13 ! runoff  
        .............
        .............
        .............
        INTEGER, PUBLIC, PARAMETER :: nDiag_ktopcl = 95 ! Level of Top of the Cloud

       ----------------- Modulo Diagnostics -----------------
      
  3.2 O segundo passo e acrescentar na variavel "avail" no nome da nova variavel
      diagnostica, lembrando que o nome nao deve passar de 40 caracteres. 
      Por exemplo:
      ----------------- Modulo Diagnostics -----------------
      avail(1:39)=(/  &
           'TIME MEAN SURFACE PRESSURE              ', &
           'TIME MEAN DIVERGENCE                    ', &
            .............
            .............
            .............
           'LEVEL OF TOP OF THE CLOUD               '/)
            
      ----------------- Modulo Diagnostics -----------------
  3.3 A variavel "lvavl" indica se a variavel diagnostica e bidimensional 
      "superficie"  ou tridimensional contendo todas as camadas verticais do
      modelo.
      Entao acrescenta 1 se a variavel diagnostica for de superficie e
      acrescenta 2 de a variavel diagnostica contiver todas as camadas 
      verticais do modelo.

 
  3.4 A variavel "nuavl" contem os codigos de conversao de unidade, caso deseja-se
      manter a unidade utilizada pelo modelo durante a computacao, entao o codido a
      ser acrescentado a esta variavel e "0" zero , Estes codigos sao indices da 
      matrix " looku " (lida do arquivo looktb), que dependendo da cobinacao 
      de indice ira fornecer os coeficientes de conversao obtidas das matrizes
      "cnfac" (lida do arquivo cnftbl) e "cnfac2" (lida do arquivo cnf2tb).

  3.4 A variavel "itavl" indica se a variavel diagnostica acrescentada e
      gaussiana ou espectral. 
      Caso deseja-se acrescentar uma varriavel diagnostica gaussiana entao
      o valor a ser acrescentado a esta variavel "itavl" e 1.
      Caso deseja-se acrescentar uma varriavel diagnostica spectral entao
      o valor a ser acrescentado a esta variavel "itavl" e 2

  3.4 A variavel "jpavl" indica em qual parte da fisica a variavel e capturada.
      1 grpcomp,physcs, 2 convection , 3 ambos, 0 neither
      Por exemplo:
      Caso deseja-se acrescentar uma varriavel diagnostica que seja computada 
      na canvection, entao o valor a ser acrescentado a esta variavel 
      "jpavl" e 2

  OBS {todas as variaveis devem ser acrescentadas no final de cada array
       descrito acima}


=========================================================================================
------------------------------------------------------
4 Module Convection:
------------------------------------------------------
  4.1 A subroutine updia e responsavel por armazenar os dados referentes
      a uma determinada variavel diagnostica.
      
      Caso deseje capturar uma determinada variavel e necessario inserir 
      a seguinte linha de codigo no local exato da captura da variavel.

      IF(dodia(nDiag_nshcrm))CALL updia(fdqn,nDiag_nshcrm,latco)
     
=========================================================================================
------------------------------------------------------
5 Tabela : desirtable        (default) 
           desirtable.pnt       
           desirtable.clm       
------------------------------------------------------
  5.1 Nas Tabelas desirtable,desirtable.pnt e desirtable.clm
      Deve acrescentar a variavel diagnostica requerida para a saida do modelo
      Po exemplo:
       ----------------- Tabela : desirtable -----------------
              40caracteres                        int  int int
        ________________________________________ ____ ____ ____
       |                                        |    |    |    |
        TOTAL PRECIPITATION                         1  121    0
        CONVECTIVE PRECIPITATION                    1  121    0
        LARGE SCALE PRECIPITATION                   1  121    0
        SNOWFALL                                    1  121    0
        RUNOFF                                      1  121    0
        INTERCEPTION LOSS                           1  170    0
        SENSIBLE HEAT FLUX FROM SURFACE             1  170    0
        LATENT HEAT FLUX FROM SURFACE               1  170    0
        SURFACE ZONAL WIND STRESS                   1  130    0
        SURFACE MERIDIONAL WIND STRESS              1  130    0
        CLOUD COVER                                 1    0    0
        DOWNWARD LONG WAVE AT BOTTOM                1  170    0
        UPWARD LONG WAVE AT BOTTOM                  1  170    0
        OUTGOING LONG WAVE AT TOP                   1  170    0
        DOWNWARD SHORT WAVE AT GROUND               1  170    0
        UPWARD SHORT WAVE AT GROUND                 1  170    0
        UPWARD SHORT WAVE AT TOP                    1  170    0
        SHORT WAVE ABSORBED AT GROUND               1  170    0
        NET LONG WAVE AT BOTTOM                     1  170    0
        GROUND/SURFACE COVER TEMPERATURE            1   40    0
        .............
        .............
        .............
        LEVEL OF TOP OF THE CLOUD                   1    0    0
      
       ----------------- Tabela : desirtable -----------------
  5.2 A primeira coluna contem o nome da variavel diagnostica, deve ter o mesmo
      nome da definicao no Modulo Diagnostic para a variavel " avail".
      A segunda coluna indica se a variavel e de superficie ou contem todos os
      niveis do modelo., deve ter o mesmo  valor da definicao no Modulo Diagnostic
      para a variavel " lvavl".
      A terceira coluna indica o codigo de conversao da variavel disgnostica, 
      deve ter o mesmo valor da definicao no Modulo Diagnostic
      para a variavel " nuavl".
      
=========================================================================================
------------------------------------------------------
6 Tabela : rfd        (default) 
           rfd.pnt       
           rfd.clm       
------------------------------------------------------
  6.1 Nas Tabelas rfd,rfd.pnt e rfd.clm
      Deve acrescentar a variavel diagnostica requerida para a saida do 
      pos-processamento
     
      Po exemplo:

       -------------------- Tabela : rfd -----------------------
              40caracteres                        int ch4 
        ________________________________________ ____ ____ ____
       |                                        |    |    |    |

        SURFACE PRESSURE                          131 PSLC  
        DIVERGENCE                                 50 DIVG
        VORTICITY                                  50 VORT
        SPECIFIC HUMIDITY                           0 UMES
        SURFACE TEMPERATURE                        40 TSFC
        SURFACE ZONAL WIND (U)                     60 UVES
        ZONAL WIND (U)                             60 UVEL
        SURFACE MERIDIONAL WIND (V)                60 VVES
        MERIDIONAL WIND (V)                        60 VVEL
        OMEGA                                     150 OMEG
        STREAM FUNCTION                            90 FCOR
        VELOCITY POTENTIAL                         90 POTV
        GEOPOTENTIAL HEIGHT                        10 ZGEO
        SEA LEVEL PRESSURE                        131 PSNM
        SURFACE ABSOLUTE TEMPERATURE               40 TEMS
        ABSOLUTE TEMPERATURE                       40 TEMP
        SURFACE RELATIVE HUMIDITY                   0 UMRS
        RELATIVE HUMIDITY                           0 UMRL
        INST. PRECIP. WATER                       110 AGPL
        TOTAL PRECIPITATION                       121 PREC
        CONVECTIVE PRECIPITATION                  121 PRCV
        OUTGOING LONG WAVE AT TOP                 170 ROLE
        CLOUD COVER                                 0 CBNV
        .............
        .............
        .............
        LEVEL OF TOP OF THE CLOUD                   0 KTOP
       -------------------- Tabela : rfd -----------------------
      
  6.2 Na tabela "rfd" estao as variaveis requerida no pos-processamento.
      O nome da variavel contida na "1 coluna" deve ser o mesmo utilizado
      pelo modelo, bem como o codigo de conversao "2 coluna" .
      O nome da variavel para a saida para o grads "3 coluna" nao ha restricao, desde
      que tenha somente 4 caracteres.
