function tab = sortTable(tab)
[tab.c,orderC] = sort(tab.c);
tab.rhohat   = tab.rhohat(orderC);
tab.kit      = tab.kit(orderC);
tab.HIV      = tab.HIV(orderC);
tab.syphilis = tab.syphilis(orderC);
tab.Ct       = tab.Ct(orderC);
tab.Ng       = tab.Ng(orderC);
end