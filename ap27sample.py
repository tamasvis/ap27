#!/usr/bin/env python3

## Demonstrate filter-bitmask construction for AP27 intervals (primes in
## a 27-element arithmetic progression).
##
## Author: Tamas Visegrady <tamas@visegrady.ch>
## SPDX-License-Identifier: BSD-3-Clause


##----------------------------------------------------------------------------
def ap27_all_terms_mod_p(a, b, p, c=223092870):
	return [ (a + b * c * i) % p  for i in range(27) ]


##----------------------------------------------------------------------------
## returns list of [ A, i ]'s where 'p' divides A+b*c*i
## i in [0, 27); A in [0, P)
##
## could be a single chained compound-compare list comprehension
## please do not comment on this
##
def ap27_composite_multiples(b, p, c=223092870):
	"A+i tuples for where AP27(A, B) has P as a small factor, or None"

	res = [ [ a, ap27_all_terms_mod_p(a, b, p, c=c) ]
	        for a in range(p) ]

			## composite if (A + B * C * i) == 0; return that i

	res = [ [ a, all.index(0) ]  for a, all in res  if (min(all) == 0) ]
	return res


##----------------------------------------------------------------------------
## returns list of mod-6*29*31 offsets, where (offset, B) may be
## starting (prime) element in an AP27 sequence
##
def ap27_mod_6_29_31(b, c=223092870):
	return [ a   for a in range(6 * 29 * 31)
	             if (((a % 6) in [ 1, 5 ]) and
	                 (min(ap27_all_terms_mod_p(a, b, 29, c=c)) >0) and
	                 (min(ap27_all_terms_mod_p(a, b, 31, c=c)) >0)) ]


##----------------------------------------------------------------------------
## pretty-print a "prime | A + (B=..) * C#23 * i" expression if known factor
##
## returns two pretty-printed terms if prime is a factor:
##     A=97*k       97|A                              OR
##     A=97*k+18    97|A + (B=1234567) * C#23 * 25
##
## returns [] if product is not known multiple of prime
##
def format_exact_division(a, b, i, p, c=223092870):
	modp = (a + b * c * i) % p
	if modp:
		return []

	afmt = f'A={p}*k' +(f'+{a % p}'  if (a % p)  else '')

	if i:
		return [ afmt, f'{p}|A + (B={b}) * C#23' +
		                         (f' * {i}'  if (i>1)  else '') ]
	return [ afmt, f'{p}|A' ]


##----------------------------------------------------------------------------
if __name__ == '__main__':
	b = 1234567

	mod1 = ap27_mod_6_29_31(b)
	print(f'## mod 6*29*31 [B={b}]; valid A(B):')
	for ai, a in enumerate(mod1):
		mod29 = ap27_all_terms_mod_p(a, b, 29)
		mod31 = ap27_all_terms_mod_p(a, b, 31)
				## AP27 mod P's; we do not log them,
				## but may be inspected here if needed
		print(f'mod 6*29*31 [{ai+1}/{len(mod1)}]', f'A={a}',
		      f'6k+{ a%6 }',
		      'mod29 all-prime:', min(mod29)>0,
		      'mod31 all-prime:', min(mod31)>0)
	print()
	assert(len(mod1) == 16)   ## by construction; sanity check


	prd = 6 * 29 * 31
	mod1offsets = sorted((prd * i) + r  for r in mod1  for i in range(32))
			## all mod-... offset in the next 32x16 offsets' list

			## show first few, original offsets
	print(f'## mod-6*29*31 offsets:', mod1offsets[:6], '...',
	      mod1offsets[-4:])
	print('##')

			## assume first step advances to first valid remainder;
			## rebase all 512 (not 511 plus a replacement tail entry),
			## so offset 0 itself is real and BIT[i] lines up with the
			## C reference's rel_offset[i] exactly, index for index
	offset1 = mod1offsets[0]
	mod1offsets = [ m - offset1  for m in mod1offsets ]

	print(f'## {len(mod1offsets)} offsets in first 32 adjacent '
		+f'mod-6*29*31 periods (first was: {offset1}; advanced) ' +
		f'[B={b}] [PRD.FULL={ 32*prd }]:')
	print(mod1offsets)
	print()

				## 53 is selected since it shows an ~even
				## mix of prohibited and possible
				## offsets (P-27 =~ P/2)
	c = 223092870
	for p in [ 37, 53, 97, 151 ]:
		print(f'## table-prime prohibited offsets({p}):')

		composite_muls = ap27_composite_multiples(b, p)

			## [ offset, AP27-index ] list for conflicts
			##
		for r in composite_muls:
			if not r:
				continue
			ai, ri = r

			aplus = f' + (B={b}) * C#23' + (f' * {ri}'  if ri>1  else '')
			##
			print(f'A={p}*k{ f"+{ai}"  if ai  else "" }    ' +
				f'{p}|A' +(aplus  if ai  else ''))
		print('##')
		print('## 512 next elems; from A === 1 mod P')

		print(f'P={p}    (B={b % p} mod P; C#23 = {c % p} mod P; '+
		      f'B * C#23 = {(b*c) %p} mod P)')

			## 512-elem bitmask for relative positions
			## show for 1 mod P as sample

		amodp, valid_bm, cmask = 1, 0, (1 << 511)
		for i, offs in enumerate(mod1offsets):
			ap = (amodp + offs) % p

			conflict = [ c  for c in composite_muls
			                if (c[0] == ap) ]

				## [ offset, AP27-index ] if conflicting
			assert(len(conflict) < 2)
				## sanity check: only single conflict (27-of-P)

			if conflict != []:
				fdiv = format_exact_division(amodp + offs, b,
							conflict[0][1], p)
				cfmt = ' ## ' + fdiv[0] +' -> ' +fdiv[1]
			else:
				valid_bm |= cmask
				cfmt = ''

			print(f'BIT[{i:3}][{offs:+7}] cand%P={ ap :<3}' + cfmt)
			cmask >>= 1

		print(f'## BITMASK[P={p}][{bin(valid_bm).count("1")}' +
			f'/{len(mod1offsets)}]=' +
			f'0x{ valid_bm :0128x}')
		print()

