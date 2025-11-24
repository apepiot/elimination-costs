function [outP_real,outP_perceived,Ph_p] = P_all(ES,kit)
%computes the asymptomatic, undiagnosed, untreated for HIV prevalence based
%on the input population of the *kit*
%% TO DO
N = sum(ES); Ph=NaN;

borneSupPrev = inf;

[HIV,syphilis,Ct,Ng,tot,~] = assigningPrevalence(ES);
if isequal(kit,{'HIV'})
    outP_real = sum(ES(HIV.asymptomaticHIV))/N;
    Ph        = outP_real;
    outP_perceived = min(f*sum(ES(HIV.asymptomaticHIV))/N,borneSupPrev);
end
if isequal(kit,{'syphilis'})
    outP_real      = sum(ES(syphilis.asymptomatic))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'Ct'})
    outP_real      = sum(ES(Ct.asymptomatic))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'Ng'})
    outP_real      = sum(ES(Ng.asymptomatic))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'HIV','syphilis'})
    asympt         = unique([HIV.asymptomaticHIV,syphilis.asymptomatic]);
    outP_real      = sum(ES(asympt))/N;
    
    Ph             = sum(ES(HIV.asymptomaticHIV))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'HIV','Ct'})
    asympt    = setdiff(unique([HIV.asymptomaticHIV,Ct.asymptomatic]),Ct.symptomatic);
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(HIV.asymptomaticHIV,Ct.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'HIV','Ng'})
    asympt = setdiff(unique([HIV.asymptomaticHIV,Ng.asymptomatic]),Ng.symptomatic);
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(HIV.asymptomaticHIV,Ng.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'syphilis','Ct'})
    asympt         = setdiff(unique([syphilis.asymptomatic,Ct.asymptomatic]),Ct.symptomatic);
    outP_real      = sum(ES(asympt))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'syphilis','Ng'})
    asympt         = setdiff(unique([syphilis.asymptomatic,Ng.asymptomatic]),Ng.symptomatic);
    outP_real      = sum(ES(asympt))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'Ct','Ng'})
    asympt         = setdiff(setdiff(unique([Ct.asymptomatic,Ng.asymptomatic]),Ct.symptomatic),Ng.symptomatic);
    outP_real      = sum(ES(asympt))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'HIV','syphilis','Ct'})
    asympt = setdiff(unique([HIV.asymptomaticHIV,syphilis.asymptomatic,Ct.asymptomatic]),Ct.symptomatic);
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(HIV.asymptomaticHIV,Ct.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'HIV','syphilis','Ng'})
    asympt = setdiff(unique([HIV.asymptomaticHIV,syphilis.asymptomatic,Ng.asymptomatic]),Ng.symptomatic);
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(HIV.asymptomaticHIV,Ng.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'HIV','Ct','Ng'})
    asympt = setdiff(setdiff(unique([HIV.asymptomaticHIV,Ct.asymptomatic,Ng.asymptomatic]),Ng.symptomatic),Ct.symptomatic);
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(setdiff(HIV.asymptomaticHIV,Ct.symptomatic),Ng.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end
if isequal(kit,{'syphilis','Ct','Ng'})
    asympt         = setdiff(setdiff(unique([syphilis.asymptomatic,Ct.asymptomatic,Ng.asymptomatic]),Ng.symptomatic),Ct.symptomatic);
    outP_real      = sum(ES(asympt))/N;
    outP_perceived = outP_real;
end
if isequal(kit,{'HIV','syphilis','Ct','Ng'})
    asympt = tot.asymptomatic;
    outP_real = sum(ES(asympt))/N;
    
    Ph             = sum(ES(setdiff(setdiff(HIV.asymptomaticHIV,Ct.symptomatic),Ng.symptomatic)))/N;
    remainedP      = outP_real - Ph;
    outP_perceived = min(f*Ph + remainedP,borneSupPrev);
end

Ph_p = f*Ph;
end


%------------------------------------------%

%     HIV1    syph1    STI1    STI2    no
%     ____    _____    ____    ____    ___
%
%     "S"     "S"      "S"     "S"       1
%     "I"     "S"      "S"     "S"       2
%     "C"     "S"      "S"     "S"       3
%     "P"     "S"      "S"     "S"       4
%     "Ip"    "S"      "S"     "S"       5
%     "Cp"    "S"      "S"     "S"       6
%     "T"     "S"      "S"     "S"       7
%     "S"     "E"      "S"     "S"       8 x
%     "I"     "E"      "S"     "S"       9
%     "C"     "E"      "S"     "S"      10
%     "P"     "E"      "S"     "S"      11 x
%     "Ip"    "E"      "S"     "S"      12
%     "Cp"    "E"      "S"     "S"      13
%     "T"     "E"      "S"     "S"      14 x
%     "S"     "I1"     "S"     "S"      15 x
%     "I"     "I1"     "S"     "S"      16
%     "C"     "I1"     "S"     "S"      17
%     "P"     "I1"     "S"     "S"      18 x
%     "Ip"    "I1"     "S"     "S"      19
%     "Cp"    "I1"     "S"     "S"      20
%     "T"     "I1"     "S"     "S"      21 x
%     "S"     "I2"     "S"     "S"      22 x
%     "I"     "I2"     "S"     "S"      23
%     "C"     "I2"     "S"     "S"      24
%     "P"     "I2"     "S"     "S"      25 x
%     "Ip"    "I2"     "S"     "S"      26
%     "Cp"    "I2"     "S"     "S"      27
%     "T"     "I2"     "S"     "S"      28 x
%     "S"     "I3"     "S"     "S"      29 x
%     "I"     "I3"     "S"     "S"      30
%     "C"     "I3"     "S"     "S"      31
%     "P"     "I3"     "S"     "S"      32 x
%     "Ip"    "I3"     "S"     "S"      33
%     "Cp"    "I3"     "S"     "S"      34
%     "T"     "I3"     "S"     "S"      35 x
%     "S"     "S"      "E"     "S"      36 x
%     "I"     "S"      "E"     "S"      37
%     "C"     "S"      "E"     "S"      38
%     "P"     "S"      "E"     "S"      39 x
%     "Ip"    "S"      "E"     "S"      40
%     "Cp"    "S"      "E"     "S"      41
%     "T"     "S"      "E"     "S"      42 x
%     "S"     "E"      "E"     "S"      43
%     "I"     "E"      "E"     "S"      44
%     "C"     "E"      "E"     "S"      45
%     "P"     "E"      "E"     "S"      46 x
%     "Ip"    "E"      "E"     "S"      47
%     "Cp"    "E"      "E"     "S"      48
%     "T"     "E"      "E"     "S"      49 x
%     "S"     "I1"     "E"     "S"      50 x
%     "I"     "I1"     "E"     "S"      51
%     "C"     "I1"     "E"     "S"      52
%     "P"     "I1"     "E"     "S"      53 x
%     "Ip"    "I1"     "E"     "S"      54
%     "Cp"    "I1"     "E"     "S"      55
%     "T"     "I1"     "E"     "S"      56 x
%     "S"     "I2"     "E"     "S"      57 x
%     "I"     "I2"     "E"     "S"      58
%     "C"     "I2"     "E"     "S"      59
%     "P"     "I2"     "E"     "S"      60 x
%     "Ip"    "I2"     "E"     "S"      61
%     "Cp"    "I2"     "E"     "S"      62
%     "T"     "I2"     "E"     "S"      63 x
%     "S"     "I3"     "E"     "S"      64
%     "I"     "I3"     "E"     "S"      65
%     "C"     "I3"     "E"     "S"      66
%     "P"     "I3"     "E"     "S"      67
%     "Ip"    "I3"     "E"     "S"      68
%     "Cp"    "I3"     "E"     "S"      69
%     "T"     "I3"     "E"     "S"      70
%     "S"     "S"      "IA"    "S"      71
%     "I"     "S"      "IA"    "S"      72
%     "C"     "S"      "IA"    "S"      73
%     "P"     "S"      "IA"    "S"      74
%     "Ip"    "S"      "IA"    "S"      75
%     "Cp"    "S"      "IA"    "S"      76
%     "T"     "S"      "IA"    "S"      77
%     "S"     "E"      "IA"    "S"      78
%     "I"     "E"      "IA"    "S"      79
%     "C"     "E"      "IA"    "S"      80
%     "P"     "E"      "IA"    "S"      81
%     "Ip"    "E"      "IA"    "S"      82
%     "Cp"    "E"      "IA"    "S"      83
%     "T"     "E"      "IA"    "S"      84
%     "S"     "I1"     "IA"    "S"      85
%     "I"     "I1"     "IA"    "S"      86
%     "C"     "I1"     "IA"    "S"      87
%     "P"     "I1"     "IA"    "S"      88
%     "Ip"    "I1"     "IA"    "S"      89
%     "Cp"    "I1"     "IA"    "S"      90
%     "T"     "I1"     "IA"    "S"      91
%     "S"     "I2"     "IA"    "S"      92
%     "I"     "I2"     "IA"    "S"      93
%     "C"     "I2"     "IA"    "S"      94
%     "P"     "I2"     "IA"    "S"      95
%     "Ip"    "I2"     "IA"    "S"      96
%     "Cp"    "I2"     "IA"    "S"      97
%     "T"     "I2"     "IA"    "S"      98
%     "S"     "I3"     "IA"    "S"      99
%     "I"     "I3"     "IA"    "S"     100
%     "C"     "I3"     "IA"    "S"     101
%     "P"     "I3"     "IA"    "S"     102
%     "Ip"    "I3"     "IA"    "S"     103
%     "Cp"    "I3"     "IA"    "S"     104
%     "T"     "I3"     "IA"    "S"     105
%     "S"     "S"      "IS"    "S"     106
%     "I"     "S"      "IS"    "S"     107
%     "C"     "S"      "IS"    "S"     108
%     "P"     "S"      "IS"    "S"     109
%     "Ip"    "S"      "IS"    "S"     110
%     "Cp"    "S"      "IS"    "S"     111
%     "T"     "S"      "IS"    "S"     112
%     "S"     "E"      "IS"    "S"     113
%     "I"     "E"      "IS"    "S"     114
%     "C"     "E"      "IS"    "S"     115
%     "P"     "E"      "IS"    "S"     116
%     "Ip"    "E"      "IS"    "S"     117
%     "Cp"    "E"      "IS"    "S"     118
%     "T"     "E"      "IS"    "S"     119
%     "S"     "I1"     "IS"    "S"     120
%     "I"     "I1"     "IS"    "S"     121
%     "C"     "I1"     "IS"    "S"     122
%     "P"     "I1"     "IS"    "S"     123
%     "Ip"    "I1"     "IS"    "S"     124
%     "Cp"    "I1"     "IS"    "S"     125
%     "T"     "I1"     "IS"    "S"     126
%     "S"     "I2"     "IS"    "S"     127
%     "I"     "I2"     "IS"    "S"     128
%     "C"     "I2"     "IS"    "S"     129
%     "P"     "I2"     "IS"    "S"     130
%     "Ip"    "I2"     "IS"    "S"     131
%     "Cp"    "I2"     "IS"    "S"     132
%     "T"     "I2"     "IS"    "S"     133
%     "S"     "I3"     "IS"    "S"     134
%     "I"     "I3"     "IS"    "S"     135
%     "C"     "I3"     "IS"    "S"     136
%     "P"     "I3"     "IS"    "S"     137
%     "Ip"    "I3"     "IS"    "S"     138
%     "Cp"    "I3"     "IS"    "S"     139
%     "T"     "I3"     "IS"    "S"     140
%     "S"     "S"      "S"     "E"     141
%     "I"     "S"      "S"     "E"     142
%     "C"     "S"      "S"     "E"     143
%     "P"     "S"      "S"     "E"     144
%     "Ip"    "S"      "S"     "E"     145
%     "Cp"    "S"      "S"     "E"     146
%     "T"     "S"      "S"     "E"     147
%     "S"     "E"      "S"     "E"     148
%     "I"     "E"      "S"     "E"     149
%     "C"     "E"      "S"     "E"     150
%     "P"     "E"      "S"     "E"     151
%     "Ip"    "E"      "S"     "E"     152
%     "Cp"    "E"      "S"     "E"     153
%     "T"     "E"      "S"     "E"     154
%     "S"     "I1"     "S"     "E"     155
%     "I"     "I1"     "S"     "E"     156
%     "C"     "I1"     "S"     "E"     157
%     "P"     "I1"     "S"     "E"     158
%     "Ip"    "I1"     "S"     "E"     159
%     "Cp"    "I1"     "S"     "E"     160
%     "T"     "I1"     "S"     "E"     161
%     "S"     "I2"     "S"     "E"     162
%     "I"     "I2"     "S"     "E"     163
%     "C"     "I2"     "S"     "E"     164
%     "P"     "I2"     "S"     "E"     165
%     "Ip"    "I2"     "S"     "E"     166
%     "Cp"    "I2"     "S"     "E"     167
%     "T"     "I2"     "S"     "E"     168
%     "S"     "I3"     "S"     "E"     169
%     "I"     "I3"     "S"     "E"     170
%     "C"     "I3"     "S"     "E"     171
%     "P"     "I3"     "S"     "E"     172
%     "Ip"    "I3"     "S"     "E"     173
%     "Cp"    "I3"     "S"     "E"     174
%     "T"     "I3"     "S"     "E"     175
%     "S"     "S"      "E"     "E"     176
%     "I"     "S"      "E"     "E"     177
%     "C"     "S"      "E"     "E"     178
%     "P"     "S"      "E"     "E"     179
%     "Ip"    "S"      "E"     "E"     180
%     "Cp"    "S"      "E"     "E"     181
%     "T"     "S"      "E"     "E"     182
%     "S"     "E"      "E"     "E"     183
%     "I"     "E"      "E"     "E"     184
%     "C"     "E"      "E"     "E"     185
%     "P"     "E"      "E"     "E"     186
%     "Ip"    "E"      "E"     "E"     187
%     "Cp"    "E"      "E"     "E"     188
%     "T"     "E"      "E"     "E"     189
%     "S"     "I1"     "E"     "E"     190
%     "I"     "I1"     "E"     "E"     191
%     "C"     "I1"     "E"     "E"     192
%     "P"     "I1"     "E"     "E"     193
%     "Ip"    "I1"     "E"     "E"     194
%     "Cp"    "I1"     "E"     "E"     195
%     "T"     "I1"     "E"     "E"     196
%     "S"     "I2"     "E"     "E"     197
%     "I"     "I2"     "E"     "E"     198
%     "C"     "I2"     "E"     "E"     199
%     "P"     "I2"     "E"     "E"     200
%     "Ip"    "I2"     "E"     "E"     201
%     "Cp"    "I2"     "E"     "E"     202
%     "T"     "I2"     "E"     "E"     203
%     "S"     "I3"     "E"     "E"     204
%     "I"     "I3"     "E"     "E"     205
%     "C"     "I3"     "E"     "E"     206
%     "P"     "I3"     "E"     "E"     207
%     "Ip"    "I3"     "E"     "E"     208
%     "Cp"    "I3"     "E"     "E"     209
%     "T"     "I3"     "E"     "E"     210
%     "S"     "S"      "IA"    "E"     211
%     "I"     "S"      "IA"    "E"     212
%     "C"     "S"      "IA"    "E"     213
%     "P"     "S"      "IA"    "E"     214
%     "Ip"    "S"      "IA"    "E"     215
%     "Cp"    "S"      "IA"    "E"     216
%     "T"     "S"      "IA"    "E"     217
%     "S"     "E"      "IA"    "E"     218
%     "I"     "E"      "IA"    "E"     219
%     "C"     "E"      "IA"    "E"     220
%     "P"     "E"      "IA"    "E"     221
%     "Ip"    "E"      "IA"    "E"     222
%     "Cp"    "E"      "IA"    "E"     223
%     "T"     "E"      "IA"    "E"     224
%     "S"     "I1"     "IA"    "E"     225
%     "I"     "I1"     "IA"    "E"     226
%     "C"     "I1"     "IA"    "E"     227
%     "P"     "I1"     "IA"    "E"     228
%     "Ip"    "I1"     "IA"    "E"     229
%     "Cp"    "I1"     "IA"    "E"     230
%     "T"     "I1"     "IA"    "E"     231
%     "S"     "I2"     "IA"    "E"     232
%     "I"     "I2"     "IA"    "E"     233
%     "C"     "I2"     "IA"    "E"     234
%     "P"     "I2"     "IA"    "E"     235
%     "Ip"    "I2"     "IA"    "E"     236
%     "Cp"    "I2"     "IA"    "E"     237
%     "T"     "I2"     "IA"    "E"     238
%     "S"     "I3"     "IA"    "E"     239
%     "I"     "I3"     "IA"    "E"     240
%     "C"     "I3"     "IA"    "E"     241
%     "P"     "I3"     "IA"    "E"     242
%     "Ip"    "I3"     "IA"    "E"     243
%     "Cp"    "I3"     "IA"    "E"     244
%     "T"     "I3"     "IA"    "E"     245
%     "S"     "S"      "IS"    "E"     246
%     "I"     "S"      "IS"    "E"     247
%     "C"     "S"      "IS"    "E"     248
%     "P"     "S"      "IS"    "E"     249
%     "Ip"    "S"      "IS"    "E"     250
%     "Cp"    "S"      "IS"    "E"     251
%     "T"     "S"      "IS"    "E"     252
%     "S"     "E"      "IS"    "E"     253
%     "I"     "E"      "IS"    "E"     254
%     "C"     "E"      "IS"    "E"     255
%     "P"     "E"      "IS"    "E"     256
%     "Ip"    "E"      "IS"    "E"     257
%     "Cp"    "E"      "IS"    "E"     258
%     "T"     "E"      "IS"    "E"     259
%     "S"     "I1"     "IS"    "E"     260
%     "I"     "I1"     "IS"    "E"     261
%     "C"     "I1"     "IS"    "E"     262
%     "P"     "I1"     "IS"    "E"     263
%     "Ip"    "I1"     "IS"    "E"     264
%     "Cp"    "I1"     "IS"    "E"     265
%     "T"     "I1"     "IS"    "E"     266
%     "S"     "I2"     "IS"    "E"     267
%     "I"     "I2"     "IS"    "E"     268
%     "C"     "I2"     "IS"    "E"     269
%     "P"     "I2"     "IS"    "E"     270
%     "Ip"    "I2"     "IS"    "E"     271
%     "Cp"    "I2"     "IS"    "E"     272
%     "T"     "I2"     "IS"    "E"     273
%     "S"     "I3"     "IS"    "E"     274
%     "I"     "I3"     "IS"    "E"     275
%     "C"     "I3"     "IS"    "E"     276
%     "P"     "I3"     "IS"    "E"     277
%     "Ip"    "I3"     "IS"    "E"     278
%     "Cp"    "I3"     "IS"    "E"     279
%     "T"     "I3"     "IS"    "E"     280
%     "S"     "S"      "S"     "IA"    281
%     "I"     "S"      "S"     "IA"    282
%     "C"     "S"      "S"     "IA"    283
%     "P"     "S"      "S"     "IA"    284
%     "Ip"    "S"      "S"     "IA"    285
%     "Cp"    "S"      "S"     "IA"    286
%     "T"     "S"      "S"     "IA"    287
%     "S"     "E"      "S"     "IA"    288
%     "I"     "E"      "S"     "IA"    289
%     "C"     "E"      "S"     "IA"    290
%     "P"     "E"      "S"     "IA"    291
%     "Ip"    "E"      "S"     "IA"    292
%     "Cp"    "E"      "S"     "IA"    293
%     "T"     "E"      "S"     "IA"    294
%     "S"     "I1"     "S"     "IA"    295
%     "I"     "I1"     "S"     "IA"    296
%     "C"     "I1"     "S"     "IA"    297
%     "P"     "I1"     "S"     "IA"    298
%     "Ip"    "I1"     "S"     "IA"    299
%     "Cp"    "I1"     "S"     "IA"    300
%     "T"     "I1"     "S"     "IA"    301
%     "S"     "I2"     "S"     "IA"    302
%     "I"     "I2"     "S"     "IA"    303
%     "C"     "I2"     "S"     "IA"    304
%     "P"     "I2"     "S"     "IA"    305
%     "Ip"    "I2"     "S"     "IA"    306
%     "Cp"    "I2"     "S"     "IA"    307
%     "T"     "I2"     "S"     "IA"    308
%     "S"     "I3"     "S"     "IA"    309
%     "I"     "I3"     "S"     "IA"    310
%     "C"     "I3"     "S"     "IA"    311
%     "P"     "I3"     "S"     "IA"    312
%     "Ip"    "I3"     "S"     "IA"    313
%     "Cp"    "I3"     "S"     "IA"    314
%     "T"     "I3"     "S"     "IA"    315
%     "S"     "S"      "E"     "IA"    316
%     "I"     "S"      "E"     "IA"    317
%     "C"     "S"      "E"     "IA"    318
%     "P"     "S"      "E"     "IA"    319
%     "Ip"    "S"      "E"     "IA"    320
%     "Cp"    "S"      "E"     "IA"    321
%     "T"     "S"      "E"     "IA"    322
%     "S"     "E"      "E"     "IA"    323
%     "I"     "E"      "E"     "IA"    324
%     "C"     "E"      "E"     "IA"    325
%     "P"     "E"      "E"     "IA"    326
%     "Ip"    "E"      "E"     "IA"    327
%     "Cp"    "E"      "E"     "IA"    328
%     "T"     "E"      "E"     "IA"    329
%     "S"     "I1"     "E"     "IA"    330
%     "I"     "I1"     "E"     "IA"    331
%     "C"     "I1"     "E"     "IA"    332
%     "P"     "I1"     "E"     "IA"    333
%     "Ip"    "I1"     "E"     "IA"    334
%     "Cp"    "I1"     "E"     "IA"    335
%     "T"     "I1"     "E"     "IA"    336
%     "S"     "I2"     "E"     "IA"    337
%     "I"     "I2"     "E"     "IA"    338
%     "C"     "I2"     "E"     "IA"    339
%     "P"     "I2"     "E"     "IA"    340
%     "Ip"    "I2"     "E"     "IA"    341
%     "Cp"    "I2"     "E"     "IA"    342
%     "T"     "I2"     "E"     "IA"    343
%     "S"     "I3"     "E"     "IA"    344
%     "I"     "I3"     "E"     "IA"    345
%     "C"     "I3"     "E"     "IA"    346
%     "P"     "I3"     "E"     "IA"    347
%     "Ip"    "I3"     "E"     "IA"    348
%     "Cp"    "I3"     "E"     "IA"    349
%     "T"     "I3"     "E"     "IA"    350
%     "S"     "S"      "IA"    "IA"    351
%     "I"     "S"      "IA"    "IA"    352
%     "C"     "S"      "IA"    "IA"    353
%     "P"     "S"      "IA"    "IA"    354
%     "Ip"    "S"      "IA"    "IA"    355
%     "Cp"    "S"      "IA"    "IA"    356
%     "T"     "S"      "IA"    "IA"    357
%     "S"     "E"      "IA"    "IA"    358
%     "I"     "E"      "IA"    "IA"    359
%     "C"     "E"      "IA"    "IA"    360
%     "P"     "E"      "IA"    "IA"    361
%     "Ip"    "E"      "IA"    "IA"    362
%     "Cp"    "E"      "IA"    "IA"    363
%     "T"     "E"      "IA"    "IA"    364
%     "S"     "I1"     "IA"    "IA"    365
%     "I"     "I1"     "IA"    "IA"    366
%     "C"     "I1"     "IA"    "IA"    367
%     "P"     "I1"     "IA"    "IA"    368
%     "Ip"    "I1"     "IA"    "IA"    369
%     "Cp"    "I1"     "IA"    "IA"    370
%     "T"     "I1"     "IA"    "IA"    371
%     "S"     "I2"     "IA"    "IA"    372
%     "I"     "I2"     "IA"    "IA"    373
%     "C"     "I2"     "IA"    "IA"    374
%     "P"     "I2"     "IA"    "IA"    375
%     "Ip"    "I2"     "IA"    "IA"    376
%     "Cp"    "I2"     "IA"    "IA"    377
%     "T"     "I2"     "IA"    "IA"    378
%     "S"     "I3"     "IA"    "IA"    379
%     "I"     "I3"     "IA"    "IA"    380
%     "C"     "I3"     "IA"    "IA"    381
%     "P"     "I3"     "IA"    "IA"    382
%     "Ip"    "I3"     "IA"    "IA"    383
%     "Cp"    "I3"     "IA"    "IA"    384
%     "T"     "I3"     "IA"    "IA"    385
%     "S"     "S"      "IS"    "IA"    386
%     "I"     "S"      "IS"    "IA"    387
%     "C"     "S"      "IS"    "IA"    388
%     "P"     "S"      "IS"    "IA"    389
%     "Ip"    "S"      "IS"    "IA"    390
%     "Cp"    "S"      "IS"    "IA"    391
%     "T"     "S"      "IS"    "IA"    392
%     "S"     "E"      "IS"    "IA"    393
%     "I"     "E"      "IS"    "IA"    394
%     "C"     "E"      "IS"    "IA"    395
%     "P"     "E"      "IS"    "IA"    396
%     "Ip"    "E"      "IS"    "IA"    397
%     "Cp"    "E"      "IS"    "IA"    398
%     "T"     "E"      "IS"    "IA"    399
%     "S"     "I1"     "IS"    "IA"    400
%     "I"     "I1"     "IS"    "IA"    401
%     "C"     "I1"     "IS"    "IA"    402
%     "P"     "I1"     "IS"    "IA"    403
%     "Ip"    "I1"     "IS"    "IA"    404
%     "Cp"    "I1"     "IS"    "IA"    405
%     "T"     "I1"     "IS"    "IA"    406
%     "S"     "I2"     "IS"    "IA"    407
%     "I"     "I2"     "IS"    "IA"    408
%     "C"     "I2"     "IS"    "IA"    409
%     "P"     "I2"     "IS"    "IA"    410
%     "Ip"    "I2"     "IS"    "IA"    411
%     "Cp"    "I2"     "IS"    "IA"    412
%     "T"     "I2"     "IS"    "IA"    413
%     "S"     "I3"     "IS"    "IA"    414
%     "I"     "I3"     "IS"    "IA"    415
%     "C"     "I3"     "IS"    "IA"    416
%     "P"     "I3"     "IS"    "IA"    417
%     "Ip"    "I3"     "IS"    "IA"    418
%     "Cp"    "I3"     "IS"    "IA"    419
%     "T"     "I3"     "IS"    "IA"    420
%     "S"     "S"      "S"     "IS"    421
%     "I"     "S"      "S"     "IS"    422
%     "C"     "S"      "S"     "IS"    423
%     "P"     "S"      "S"     "IS"    424
%     "Ip"    "S"      "S"     "IS"    425
%     "Cp"    "S"      "S"     "IS"    426
%     "T"     "S"      "S"     "IS"    427
%     "S"     "E"      "S"     "IS"    428
%     "I"     "E"      "S"     "IS"    429
%     "C"     "E"      "S"     "IS"    430
%     "P"     "E"      "S"     "IS"    431
%     "Ip"    "E"      "S"     "IS"    432
%     "Cp"    "E"      "S"     "IS"    433
%     "T"     "E"      "S"     "IS"    434
%     "S"     "I1"     "S"     "IS"    435
%     "I"     "I1"     "S"     "IS"    436
%     "C"     "I1"     "S"     "IS"    437
%     "P"     "I1"     "S"     "IS"    438
%     "Ip"    "I1"     "S"     "IS"    439
%     "Cp"    "I1"     "S"     "IS"    440
%     "T"     "I1"     "S"     "IS"    441
%     "S"     "I2"     "S"     "IS"    442
%     "I"     "I2"     "S"     "IS"    443
%     "C"     "I2"     "S"     "IS"    444
%     "P"     "I2"     "S"     "IS"    445
%     "Ip"    "I2"     "S"     "IS"    446
%     "Cp"    "I2"     "S"     "IS"    447
%     "T"     "I2"     "S"     "IS"    448
%     "S"     "I3"     "S"     "IS"    449
%     "I"     "I3"     "S"     "IS"    450
%     "C"     "I3"     "S"     "IS"    451
%     "P"     "I3"     "S"     "IS"    452
%     "Ip"    "I3"     "S"     "IS"    453
%     "Cp"    "I3"     "S"     "IS"    454
%     "T"     "I3"     "S"     "IS"    455
%     "S"     "S"      "E"     "IS"    456
%     "I"     "S"      "E"     "IS"    457
%     "C"     "S"      "E"     "IS"    458
%     "P"     "S"      "E"     "IS"    459
%     "Ip"    "S"      "E"     "IS"    460
%     "Cp"    "S"      "E"     "IS"    461
%     "T"     "S"      "E"     "IS"    462
%     "S"     "E"      "E"     "IS"    463
%     "I"     "E"      "E"     "IS"    464
%     "C"     "E"      "E"     "IS"    465
%     "P"     "E"      "E"     "IS"    466
%     "Ip"    "E"      "E"     "IS"    467
%     "Cp"    "E"      "E"     "IS"    468
%     "T"     "E"      "E"     "IS"    469
%     "S"     "I1"     "E"     "IS"    470
%     "I"     "I1"     "E"     "IS"    471
%     "C"     "I1"     "E"     "IS"    472
%     "P"     "I1"     "E"     "IS"    473
%     "Ip"    "I1"     "E"     "IS"    474
%     "Cp"    "I1"     "E"     "IS"    475
%     "T"     "I1"     "E"     "IS"    476
%     "S"     "I2"     "E"     "IS"    477
%     "I"     "I2"     "E"     "IS"    478
%     "C"     "I2"     "E"     "IS"    479
%     "P"     "I2"     "E"     "IS"    480
%     "Ip"    "I2"     "E"     "IS"    481
%     "Cp"    "I2"     "E"     "IS"    482
%     "T"     "I2"     "E"     "IS"    483
%     "S"     "I3"     "E"     "IS"    484
%     "I"     "I3"     "E"     "IS"    485
%     "C"     "I3"     "E"     "IS"    486
%     "P"     "I3"     "E"     "IS"    487
%     "Ip"    "I3"     "E"     "IS"    488
%     "Cp"    "I3"     "E"     "IS"    489
%     "T"     "I3"     "E"     "IS"    490
%     "S"     "S"      "IA"    "IS"    491
%     "I"     "S"      "IA"    "IS"    492
%     "C"     "S"      "IA"    "IS"    493
%     "P"     "S"      "IA"    "IS"    494
%     "Ip"    "S"      "IA"    "IS"    495
%     "Cp"    "S"      "IA"    "IS"    496
%     "T"     "S"      "IA"    "IS"    497
%     "S"     "E"      "IA"    "IS"    498
%     "I"     "E"      "IA"    "IS"    499
%     "C"     "E"      "IA"    "IS"    500
%     "P"     "E"      "IA"    "IS"    501
%     "Ip"    "E"      "IA"    "IS"    502
%     "Cp"    "E"      "IA"    "IS"    503
%     "T"     "E"      "IA"    "IS"    504
%     "S"     "I1"     "IA"    "IS"    505
%     "I"     "I1"     "IA"    "IS"    506
%     "C"     "I1"     "IA"    "IS"    507
%     "P"     "I1"     "IA"    "IS"    508
%     "Ip"    "I1"     "IA"    "IS"    509
%     "Cp"    "I1"     "IA"    "IS"    510
%     "T"     "I1"     "IA"    "IS"    511
%     "S"     "I2"     "IA"    "IS"    512
%     "I"     "I2"     "IA"    "IS"    513
%     "C"     "I2"     "IA"    "IS"    514
%     "P"     "I2"     "IA"    "IS"    515
%     "Ip"    "I2"     "IA"    "IS"    516
%     "Cp"    "I2"     "IA"    "IS"    517
%     "T"     "I2"     "IA"    "IS"    518
%     "S"     "I3"     "IA"    "IS"    519
%     "I"     "I3"     "IA"    "IS"    520
%     "C"     "I3"     "IA"    "IS"    521
%     "P"     "I3"     "IA"    "IS"    522
%     "Ip"    "I3"     "IA"    "IS"    523
%     "Cp"    "I3"     "IA"    "IS"    524
%     "T"     "I3"     "IA"    "IS"    525
%     "S"     "S"      "IS"    "IS"    526
%     "I"     "S"      "IS"    "IS"    527
%     "C"     "S"      "IS"    "IS"    528
%     "P"     "S"      "IS"    "IS"    529
%     "Ip"    "S"      "IS"    "IS"    530
%     "Cp"    "S"      "IS"    "IS"    531
%     "T"     "S"      "IS"    "IS"    532
%     "S"     "E"      "IS"    "IS"    533
%     "I"     "E"      "IS"    "IS"    534
%     "C"     "E"      "IS"    "IS"    535
%     "P"     "E"      "IS"    "IS"    536
%     "Ip"    "E"      "IS"    "IS"    537
%     "Cp"    "E"      "IS"    "IS"    538
%     "T"     "E"      "IS"    "IS"    539
%     "S"     "I1"     "IS"    "IS"    540
%     "I"     "I1"     "IS"    "IS"    541
%     "C"     "I1"     "IS"    "IS"    542
%     "P"     "I1"     "IS"    "IS"    543
%     "Ip"    "I1"     "IS"    "IS"    544
%     "Cp"    "I1"     "IS"    "IS"    545
%     "T"     "I1"     "IS"    "IS"    546
%     "S"     "I2"     "IS"    "IS"    547
%     "I"     "I2"     "IS"    "IS"    548
%     "C"     "I2"     "IS"    "IS"    549
%     "P"     "I2"     "IS"    "IS"    550
%     "Ip"    "I2"     "IS"    "IS"    551
%     "Cp"    "I2"     "IS"    "IS"    552
%     "T"     "I2"     "IS"    "IS"    553
%     "S"     "I3"     "IS"    "IS"    554
%     "I"     "I3"     "IS"    "IS"    555
%     "C"     "I3"     "IS"    "IS"    556
%     "P"     "I3"     "IS"    "IS"    557
%     "Ip"    "I3"     "IS"    "IS"    558
%     "Cp"    "I3"     "IS"    "IS"    559
%     "T"     "I3"     "IS"    "IS"    560