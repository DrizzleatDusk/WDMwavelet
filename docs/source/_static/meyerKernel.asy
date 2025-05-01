size(10cm);
defaultpen(fontsize(10pt));

int Nt = 5;
real A = 1.5;
real B = 2.;
real H = .8 * Nt;
real marg = .04 * Nt;

Label lx = Label("$k$", EndPoint, align=S);
Label ly = Label("$\tilde{\phi}$", EndPoint, align=W);

label("$O$", (0, 0), align=SW);
draw((-Nt-1, 0)--(6, 0));
draw((6.2, 0)--(2.55*Nt, 0), arrow=ArcArrow(SimpleHead), L=lx);
draw((0, -Nt/4)--(0, Nt), arrow=ArcArrow(SimpleHead), L=ly);

path ker_r = (0, H)--(A, H){(1, 0)}..(A+B/4, .97*H)..(A+B/2, .75*H)..(A+3*B/4, .15*H)..{(1, 0)}(A+B, 0)--(Nt, 0);
path ker_l = reflect((0, 0), (0, 1)) * ker_r;
//path ker_ls = shift() * ker_l;
draw(ker_r, blue);
draw((A, 0)--(A, H), dashed);
draw((2.4*Nt, H)--(2.4*Nt, 0), dashed);
draw(ker_l, blue+opacity(.4));
dot((0, H), blue);
dot((1, H), blue);
dot((2, .965*H), blue);
dot((3, .15*H), blue);
dot((4, 0), blue);
dot((5, 0), blue);


Label lA = Label("$A$", MidPoint, align=N);
Label lB = Label("$B$", MidPoint, align=N);
draw((0, marg)--(A, marg), arrow=Arrows(6), bar=Bars, L=lA);
draw((A, marg)--(A+B, marg), arrow=Arrows(6), bar=Bars, L=lB);
draw((A+B, marg)--(2*A+B, marg), arrow=Arrows(6), bar=Bars, L=lA);

path ker_ls = reverse(shift((2.4*Nt, 0)) * ker_l);
draw(ker_ls, blue, arrow=ArcArrow(SimpleHead));

Label lNt = Label("$N_t$", align=S);
Label lN = Label("$N$", align=S);
dot((Nt, 0), L=lNt);
dot((5, 0), blue);
dot((2.4*Nt, 0), L=lN);

draw((6+marg, 1.8*marg)--(6-marg, -1.8*marg));
draw((6.2+marg, 1.8*marg)--(6.2-marg, -1.8*marg));