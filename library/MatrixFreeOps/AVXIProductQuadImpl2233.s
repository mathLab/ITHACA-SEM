00000000000b92d0 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE>:
   b92d0:   4c 8d 54 24 08          lea    0x8(%rsp),%r10
   b92d5:   48 83 e4 e0             and    $0xffffffffffffffe0,%rsp
   b92d9:   41 ff 72 f8             pushq  -0x8(%r10)
   b92dd:   55                      push   %rbp
   b92de:   48 89 e5                mov    %rsp,%rbp
   b92e1:   41 57                   push   %r15
   b92e3:   41 56                   push   %r14
   b92e5:   41 55                   push   %r13
   b92e7:   41 54                   push   %r12
   b92e9:   49 89 fc                mov    %rdi,%r12
   b92ec:   41 52                   push   %r10
   b92ee:   53                      push   %rbx
   b92ef:   48 83 c4 80             add    $0xffffffffffffff80,%rsp
   b92f3:   48 8b 4e 20             mov    0x20(%rsi),%rcx
   b92f7:   64 48 8b 04 25 28 00    mov    %fs:0x28,%rax
   b92fe:   00 00
   b9300:   48 89 45 c8             mov    %rax,-0x38(%rbp)
   b9304:   31 c0                   xor    %eax,%eax
   b9306:   48 8b 46 10             mov    0x10(%rsi),%rax
   b930a:   be 09 00 00 00          mov    $0x9,%esi
   b930f:   48 8d 1c c8             lea    (%rax,%rcx,8),%rbx
   b9313:   48 8b 42 10             mov    0x10(%rdx),%rax
   b9317:   48 8b 4a 20             mov    0x20(%rdx),%rcx
   b931b:   4c 8d 2c c8             lea    (%rax,%rcx,8),%r13
   b931f:   8b 87 58 02 00 00       mov    0x258(%rdi),%eax
   b9325:   48 8d 7d 90             lea    -0x70(%rbp),%rdi
   b9329:   44 8d 34 85 00 00 00    lea    0x0(,%rax,4),%r14d
   b9330:   00
   b9331:   e8 fa c6 fe ff          callq  a5a30 <_ZNSt6vectorIN6Nektar3AVX7VecDataIdLi4EEEN5boost9alignment17aligned_allocatorIS3_Lm32EEEEC2EmRKS7_.constprop.404>
   b9336:   49 63 b4 24 58 02 00    movslq 0x258(%r12),%rsi
   b933d:   00
   b933e:   48 8d 45 b0             lea    -0x50(%rbp),%rax
   b9342:   48 89 c7                mov    %rax,%rdi
   b9345:   48 89 85 68 ff ff ff    mov    %rax,-0x98(%rbp)
   b934c:   e8 df c6 fe ff          callq  a5a30 <_ZNSt6vectorIN6Nektar3AVX7VecDataIdLi4EEEN5boost9alignment17aligned_allocatorIS3_Lm32EEEEC2EmRKS7_.constprop.404>
   b9351:   41 8b 44 24 30          mov    0x30(%r12),%eax
   b9356:   85 c0                   test   %eax,%eax
   b9358:   0f 8e a6 03 00 00       jle    b9704 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x434>
   b935e:   4d 63 f6                movslq %r14d,%r14
   b9361:   c5 e1 57 db             vxorpd %xmm3,%xmm3,%xmm3
   b9365:   4a 8d 04 f5 00 00 00    lea    0x0(,%r14,8),%rax
   b936c:   00
   b936d:   45 31 ff                xor    %r15d,%r15d
   b9370:   45 31 f6                xor    %r14d,%r14d
   b9373:   48 89 85 60 ff ff ff    mov    %rax,-0xa0(%rbp)
   b937a:   66 0f 1f 44 00 00       nopw   0x0(%rax,%rax,1)
   b9380:   c5 f9 10 13             vmovupd (%rbx),%xmm2
   b9384:   48 8b 45 90             mov    -0x70(%rbp),%rax
   b9388:   4d 8b 9c 24 50 01 00    mov    0x150(%r12),%r11
   b938f:   00
   b9390:   c5 f9 10 6b 48          vmovupd 0x48(%rbx),%xmm5
   b9395:   c4 e3 6d 18 93 90 00    vinsertf128 $0x1,0x90(%rbx),%ymm2,%ymm2
   b939c:   00 00 01
   b939f:   c5 f9 10 4b 10          vmovupd 0x10(%rbx),%xmm1
   b93a4:   4d 01 fb                add    %r15,%r11
   b93a7:   c4 e3 55 18 ab d8 00    vinsertf128 $0x1,0xd8(%rbx),%ymm5,%ymm5
   b93ae:   00 00 01
   b93b1:   c5 f9 10 63 58          vmovupd 0x58(%rbx),%xmm4
   b93b6:   c4 e3 75 18 8b a0 00    vinsertf128 $0x1,0xa0(%rbx),%ymm1,%ymm1
   b93bd:   00 00 01
   b93c0:   c4 e3 5d 18 a3 e8 00    vinsertf128 $0x1,0xe8(%rbx),%ymm4,%ymm4
   b93c7:   00 00 01
   b93ca:   c5 ed 14 c5             vunpcklpd %ymm5,%ymm2,%ymm0
   b93ce:   c5 ed 15 d5             vunpckhpd %ymm5,%ymm2,%ymm2
   b93d2:   c5 fd 29 00             vmovapd %ymm0,(%rax)
   b93d6:   c5 fd 29 50 20          vmovapd %ymm2,0x20(%rax)
   b93db:   c5 f5 14 d4             vunpcklpd %ymm4,%ymm1,%ymm2
   b93df:   c5 f5 15 cc             vunpckhpd %ymm4,%ymm1,%ymm1
   b93e3:   c5 fd 29 50 40          vmovapd %ymm2,0x40(%rax)
   b93e8:   c5 fd 29 48 60          vmovapd %ymm1,0x60(%rax)
   b93ed:   c5 f9 10 53 20          vmovupd 0x20(%rbx),%xmm2
   b93f2:   c5 f9 10 6b 68          vmovupd 0x68(%rbx),%xmm5
   b93f7:   c4 e3 6d 18 93 b0 00    vinsertf128 $0x1,0xb0(%rbx),%ymm2,%ymm2
   b93fe:   00 00 01
   b9401:   c5 f9 10 4b 30          vmovupd 0x30(%rbx),%xmm1
   b9406:   c4 e3 55 18 ab f8 00    vinsertf128 $0x1,0xf8(%rbx),%ymm5,%ymm5
   b940d:   00 00 01
   b9410:   c5 f9 10 63 78          vmovupd 0x78(%rbx),%xmm4
   b9415:   c4 e3 75 18 8b c0 00    vinsertf128 $0x1,0xc0(%rbx),%ymm1,%ymm1
   b941c:   00 00 01
   b941f:   c4 e3 5d 18 a3 08 01    vinsertf128 $0x1,0x108(%rbx),%ymm4,%ymm4
   b9426:   00 00 01
   b9429:   c5 ed 14 f5             vunpcklpd %ymm5,%ymm2,%ymm6
   b942d:   c5 ed 15 d5             vunpckhpd %ymm5,%ymm2,%ymm2
   b9431:   c5 fd 29 b0 80 00 00    vmovapd %ymm6,0x80(%rax)
   b9438:   00
   b9439:   c5 fd 29 90 a0 00 00    vmovapd %ymm2,0xa0(%rax)
   b9440:   00
   b9441:   c5 f5 14 d4             vunpcklpd %ymm4,%ymm1,%ymm2
   b9445:   c5 f5 15 cc             vunpckhpd %ymm4,%ymm1,%ymm1
   b9449:   c5 fd 29 90 c0 00 00    vmovapd %ymm2,0xc0(%rax)
   b9450:   00
   b9451:   c5 fd 29 88 e0 00 00    vmovapd %ymm1,0xe0(%rax)
   b9458:   00
   b9459:   c5 fb 10 4b 40          vmovsd 0x40(%rbx),%xmm1
   b945e:   c5 fb 11 88 00 01 00    vmovsd %xmm1,0x100(%rax)
   b9465:   00
   b9466:   c5 fb 10 8b 88 00 00    vmovsd 0x88(%rbx),%xmm1
   b946d:   00
   b946e:   c5 fb 11 88 08 01 00    vmovsd %xmm1,0x108(%rax)
   b9475:   00
   b9476:   c5 fb 10 8b d0 00 00    vmovsd 0xd0(%rbx),%xmm1
   b947d:   00
   b947e:   c5 fb 11 88 10 01 00    vmovsd %xmm1,0x110(%rax)
   b9485:   00
   b9486:   c5 fb 10 8b 18 01 00    vmovsd 0x118(%rbx),%xmm1
   b948d:   00
   b948e:   49 8b 7c 24 48          mov    0x48(%r12),%rdi
   b9493:   49 8b b4 24 08 01 00    mov    0x108(%r12),%rsi
   b949a:   00
   b949b:   c5 fb 11 88 18 01 00    vmovsd %xmm1,0x118(%rax)
   b94a2:   00
   b94a3:   49 8b 4c 24 60          mov    0x60(%r12),%rcx
   b94a8:   49 8b 94 24 20 01 00    mov    0x120(%r12),%rdx
   b94af:   00
   b94b0:   4c 8b 55 b0             mov    -0x50(%rbp),%r10
   b94b4:   c5 fd 28 2f             vmovapd (%rdi),%ymm5
   b94b8:   c5 fd 28 57 20          vmovapd 0x20(%rdi),%ymm2
   b94bd:   c5 fd 59 cd             vmulpd %ymm5,%ymm0,%ymm1
   b94c1:   c4 c1 7d 28 23          vmovapd (%r11),%ymm4
   b94c6:   c5 7d 28 0e             vmovapd (%rsi),%ymm9
   b94ca:   c5 ed 59 40 20          vmulpd 0x20(%rax),%ymm2,%ymm0
   b94cf:   c5 7d 28 46 20          vmovapd 0x20(%rsi),%ymm8
   b94d4:   c5 fd 28 7e 40          vmovapd 0x40(%rsi),%ymm7
   b94d9:   c5 f5 59 cc             vmulpd %ymm4,%ymm1,%ymm1
   b94dd:   c5 55 59 50 60          vmulpd 0x60(%rax),%ymm5,%ymm10
   b94e2:   c5 fd 59 c4             vmulpd %ymm4,%ymm0,%ymm0
   b94e6:   c4 c2 e5 98 c9          vfmadd132pd %ymm9,%ymm3,%ymm1
   b94eb:   c5 d5 59 a8 c0 00 00    vmulpd 0xc0(%rax),%ymm5,%ymm5
   b94f2:   00
   b94f3:   c4 c2 fd b8 c8          vfmadd231pd %ymm8,%ymm0,%ymm1
   b94f8:   c5 fd 28 47 40          vmovapd 0x40(%rdi),%ymm0
   b94fd:   c5 2d 59 d4             vmulpd %ymm4,%ymm10,%ymm10
   b9501:   c5 fd 59 70 40          vmulpd 0x40(%rax),%ymm0,%ymm6
   b9506:   c5 d5 59 ec             vmulpd %ymm4,%ymm5,%ymm5
   b950a:   c4 42 e5 98 d1          vfmadd132pd %ymm9,%ymm3,%ymm10
   b950f:   c5 cd 59 f4             vmulpd %ymm4,%ymm6,%ymm6
   b9513:   c4 c2 e5 98 e9          vfmadd132pd %ymm9,%ymm3,%ymm5
   b9518:   c4 e2 cd b8 cf          vfmadd231pd %ymm7,%ymm6,%ymm1
   b951d:   c5 ed 59 b0 80 00 00    vmulpd 0x80(%rax),%ymm2,%ymm6
   b9524:   00
   b9525:   c5 ed 59 90 e0 00 00    vmulpd 0xe0(%rax),%ymm2,%ymm2
   b952c:   00
   b952d:   c5 cd 59 f4             vmulpd %ymm4,%ymm6,%ymm6
   b9531:   c5 ed 59 d4             vmulpd %ymm4,%ymm2,%ymm2
   b9535:   c4 42 cd b8 d0          vfmadd231pd %ymm8,%ymm6,%ymm10
   b953a:   c5 fd 59 b0 a0 00 00    vmulpd 0xa0(%rax),%ymm0,%ymm6
   b9541:   00
   b9542:   c5 fd 59 80 00 01 00    vmulpd 0x100(%rax),%ymm0,%ymm0
   b9549:   00
   b954a:   c4 c2 d5 98 d0          vfmadd132pd %ymm8,%ymm5,%ymm2
   b954f:   c5 cd 59 f4             vmulpd %ymm4,%ymm6,%ymm6
   b9553:   c5 fd 59 c4             vmulpd %ymm4,%ymm0,%ymm0
   b9557:   c5 f5 59 21             vmulpd (%rcx),%ymm1,%ymm4
   b955b:   c4 e2 ad 98 f7          vfmadd132pd %ymm7,%ymm10,%ymm6
   b9560:   c4 e2 ed 98 c7          vfmadd132pd %ymm7,%ymm2,%ymm0
   b9565:   c5 cd 59 51 20          vmulpd 0x20(%rcx),%ymm6,%ymm2
   b956a:   c4 e2 e5 98 22          vfmadd132pd (%rdx),%ymm3,%ymm4
   b956f:   c4 e2 ed b8 62 20       vfmadd231pd 0x20(%rdx),%ymm2,%ymm4
   b9575:   c5 fd 59 51 40          vmulpd 0x40(%rcx),%ymm0,%ymm2
   b957a:   c4 e2 dd 98 52 40       vfmadd132pd 0x40(%rdx),%ymm4,%ymm2
   b9580:   c4 c1 7d 29 12          vmovapd %ymm2,(%r10)
   b9585:   c5 f5 59 49 60          vmulpd 0x60(%rcx),%ymm1,%ymm1
   b958a:   c5 cd 59 b1 80 00 00    vmulpd 0x80(%rcx),%ymm6,%ymm6
   b9591:   00
   b9592:   c4 e2 e5 98 0a          vfmadd132pd (%rdx),%ymm3,%ymm1
   b9597:   c4 e2 f5 98 72 20       vfmadd132pd 0x20(%rdx),%ymm1,%ymm6
   b959d:   c5 fd 59 81 a0 00 00    vmulpd 0xa0(%rcx),%ymm0,%ymm0
   b95a4:   00
   b95a5:   c5 fd 29 9d 70 ff ff    vmovapd %ymm3,-0x90(%rbp)
   b95ac:   ff
   b95ad:   c4 e2 cd 98 42 40       vfmadd132pd 0x40(%rdx),%ymm6,%ymm0
   b95b3:   c4 c1 7d 29 42 40       vmovapd %ymm0,0x40(%r10)
   b95b9:   c5 fd 28 6f 60          vmovapd 0x60(%rdi),%ymm5
   b95be:   c5 fd 28 a7 80 00 00    vmovapd 0x80(%rdi),%ymm4
   b95c5:   00
   b95c6:   c5 d5 59 08             vmulpd (%rax),%ymm5,%ymm1
   b95ca:   c4 c1 7d 28 33          vmovapd (%r11),%ymm6
   b95cf:   c5 dd 59 50 20          vmulpd 0x20(%rax),%ymm4,%ymm2
   b95d4:   c5 7d 28 16             vmovapd (%rsi),%ymm10
   b95d8:   c5 7d 28 4e 20          vmovapd 0x20(%rsi),%ymm9
   b95dd:   c5 fd 28 87 a0 00 00    vmovapd 0xa0(%rdi),%ymm0
   b95e4:   00
   b95e5:   c5 fd 28 7e 40          vmovapd 0x40(%rsi),%ymm7
   b95ea:   c5 f5 59 ce             vmulpd %ymm6,%ymm1,%ymm1
   b95ee:   c5 ed 59 d6             vmulpd %ymm6,%ymm2,%ymm2
   b95f2:   c5 55 59 40 60          vmulpd 0x60(%rax),%ymm5,%ymm8
   b95f7:   c4 c2 e5 98 ca          vfmadd132pd %ymm10,%ymm3,%ymm1
   b95fc:   c5 d5 59 a8 c0 00 00    vmulpd 0xc0(%rax),%ymm5,%ymm5
   b9603:   00
   b9604:   c4 c2 ed b8 c9          vfmadd231pd %ymm9,%ymm2,%ymm1
   b9609:   c5 fd 59 50 40          vmulpd 0x40(%rax),%ymm0,%ymm2
   b960e:   c5 3d 59 c6             vmulpd %ymm6,%ymm8,%ymm8
   b9612:   c5 d5 59 ee             vmulpd %ymm6,%ymm5,%ymm5
   b9616:   c5 ed 59 d6             vmulpd %ymm6,%ymm2,%ymm2
   b961a:   c4 42 e5 98 c2          vfmadd132pd %ymm10,%ymm3,%ymm8
   b961f:   c4 c2 e5 98 ea          vfmadd132pd %ymm10,%ymm3,%ymm5
   b9624:   c4 e2 f5 98 d7          vfmadd132pd %ymm7,%ymm1,%ymm2
   b9629:   c5 dd 59 88 80 00 00    vmulpd 0x80(%rax),%ymm4,%ymm1
   b9630:   00
   b9631:   c5 dd 59 a0 e0 00 00    vmulpd 0xe0(%rax),%ymm4,%ymm4
   b9638:   00
   b9639:   c5 f5 59 ce             vmulpd %ymm6,%ymm1,%ymm1
   b963d:   c5 dd 59 e6             vmulpd %ymm6,%ymm4,%ymm4
   b9641:   c4 42 f5 b8 c1          vfmadd231pd %ymm9,%ymm1,%ymm8
   b9646:   c5 fd 59 88 a0 00 00    vmulpd 0xa0(%rax),%ymm0,%ymm1
   b964d:   00
   b964e:   c5 fd 59 80 00 01 00    vmulpd 0x100(%rax),%ymm0,%ymm0
   b9655:   00
   b9656:   c4 c2 d5 98 e1          vfmadd132pd %ymm9,%ymm5,%ymm4
   b965b:   c5 ed 59 29             vmulpd (%rcx),%ymm2,%ymm5
   b965f:   c5 f5 59 ce             vmulpd %ymm6,%ymm1,%ymm1
   b9663:   c5 fd 59 c6             vmulpd %ymm6,%ymm0,%ymm0
   b9667:   c4 e2 e5 98 2a          vfmadd132pd (%rdx),%ymm3,%ymm5
   b966c:   c4 e2 bd 98 cf          vfmadd132pd %ymm7,%ymm8,%ymm1
   b9671:   c4 e2 dd 98 c7          vfmadd132pd %ymm7,%ymm4,%ymm0
   b9676:   c5 f5 59 61 20          vmulpd 0x20(%rcx),%ymm1,%ymm4
   b967b:   c4 e2 dd b8 6a 20       vfmadd231pd 0x20(%rdx),%ymm4,%ymm5
   b9681:   c5 fd 59 61 40          vmulpd 0x40(%rcx),%ymm0,%ymm4
   b9686:   c4 e2 d5 98 62 40       vfmadd132pd 0x40(%rdx),%ymm5,%ymm4
   b968c:   c4 c1 7d 29 62 20       vmovapd %ymm4,0x20(%r10)
   b9692:   c5 ed 59 51 60          vmulpd 0x60(%rcx),%ymm2,%ymm2
   b9697:   c5 f5 59 89 80 00 00    vmulpd 0x80(%rcx),%ymm1,%ymm1
   b969e:   00
   b969f:   c4 e2 e5 98 12          vfmadd132pd (%rdx),%ymm3,%ymm2
   b96a4:   c4 e2 ed 98 4a 20       vfmadd132pd 0x20(%rdx),%ymm2,%ymm1
   b96aa:   c5 fd 59 81 a0 00 00    vmulpd 0xa0(%rcx),%ymm0,%ymm0
   b96b1:   00
   b96b2:   49 63 b4 24 58 02 00    movslq 0x258(%r12),%rsi
   b96b9:   00
   b96ba:   48 8b bd 68 ff ff ff    mov    -0x98(%rbp),%rdi
   b96c1:   c4 e2 f5 98 42 40       vfmadd132pd 0x40(%rdx),%ymm1,%ymm0
   b96c7:   4c 89 ea                mov    %r13,%rdx
   b96ca:   c4 c1 7d 29 42 60       vmovapd %ymm0,0x60(%r10)
   b96d0:   c5 f8 77                vzeroupper
   b96d3:   e8 98 db f6 ff          callq  27270 <_ZN6Nektar3AVX7VecDataIdLi4EE18deinterleave_storeERKSt6vectorIS2_N5boost9alignment17aligned_allocatorIS2_Lm32EEEEmPd@plt>
   b96d8:   48 81 c3 20 01 00 00    add    $0x120,%rbx
   b96df:   4c 03 ad 60 ff ff ff    add    -0xa0(%rbp),%r13
   b96e6:   41 83 c6 01             add    $0x1,%r14d
   b96ea:   49 83 c7 20             add    $0x20,%r15
   b96ee:   45 39 74 24 30          cmp    %r14d,0x30(%r12)
   b96f3:   c5 fd 28 9d 70 ff ff    vmovapd -0x90(%rbp),%ymm3
   b96fa:   ff
   b96fb:   0f 8f 7f fc ff ff       jg     b9380 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0xb0>
   b9701:   c5 f8 77                vzeroupper
   b9704:   48 8b 7d b0             mov    -0x50(%rbp),%rdi
   b9708:   48 85 ff                test   %rdi,%rdi
   b970b:   74 05                   je     b9712 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x442>
   b970d:   e8 de d3 f6 ff          callq  26af0 <free@plt>
   b9712:   48 8b 7d 90             mov    -0x70(%rbp),%rdi
   b9716:   48 85 ff                test   %rdi,%rdi
   b9719:   74 05                   je     b9720 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x450>
   b971b:   e8 d0 d3 f6 ff          callq  26af0 <free@plt>
   b9720:   48 8b 45 c8             mov    -0x38(%rbp),%rax
   b9724:   64 48 33 04 25 28 00    xor    %fs:0x28,%rax
   b972b:   00 00
   b972d:   75 31                   jne    b9760 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x490>
   b972f:   48 83 ec 80             sub    $0xffffffffffffff80,%rsp
   b9733:   5b                      pop    %rbx
   b9734:   41 5a                   pop    %r10
   b9736:   41 5c                   pop    %r12
   b9738:   41 5d                   pop    %r13
   b973a:   41 5e                   pop    %r14
   b973c:   41 5f                   pop    %r15
   b973e:   5d                      pop    %rbp
   b973f:   49 8d 62 f8             lea    -0x8(%r10),%rsp
   b9743:   c3                      retq
   b9744:   48 8b 7d 90             mov    -0x70(%rbp),%rdi
   b9748:   48 89 c3                mov    %rax,%rbx
   b974b:   48 85 ff                test   %rdi,%rdi
   b974e:   74 15                   je     b9765 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x495>
   b9750:   c5 f8 77                vzeroupper
   b9753:   e8 98 d3 f6 ff          callq  26af0 <free@plt>
   b9758:   48 89 df                mov    %rbx,%rdi
   b975b:   e8 e0 da f6 ff          callq  27240 <_Unwind_Resume@plt>
   b9760:   e8 0b d2 f6 ff          callq  26970 <__stack_chk_fail@plt>
   b9765:   c5 f8 77                vzeroupper
   b9768:   eb ee                   jmp    b9758 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb0EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x488>
   b976a:   66 0f 1f 44 00 00       nopw   0x0(%rax,%rax,1)

00000000000aab10 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE>:
   aab10:   4c 8d 54 24 08          lea    0x8(%rsp),%r10
   aab15:   48 83 e4 e0             and    $0xffffffffffffffe0,%rsp
   aab19:   41 ff 72 f8             pushq  -0x8(%r10)
   aab1d:   55                      push   %rbp
   aab1e:   48 89 e5                mov    %rsp,%rbp
   aab21:   41 57                   push   %r15
   aab23:   41 56                   push   %r14
   aab25:   41 55                   push   %r13
   aab27:   41 54                   push   %r12
   aab29:   41 52                   push   %r10
   aab2b:   53                      push   %rbx
   aab2c:   48 89 fb                mov    %rdi,%rbx
   aab2f:   48 83 c4 80             add    $0xffffffffffffff80,%rsp
   aab33:   48 8b 4e 20             mov    0x20(%rsi),%rcx
   aab37:   64 48 8b 04 25 28 00    mov    %fs:0x28,%rax
   aab3e:   00 00
   aab40:   48 89 45 c8             mov    %rax,-0x38(%rbp)
   aab44:   31 c0                   xor    %eax,%eax
   aab46:   48 8b 46 10             mov    0x10(%rsi),%rax
   aab4a:   be 09 00 00 00          mov    $0x9,%esi
   aab4f:   4c 8d 3c c8             lea    (%rax,%rcx,8),%r15
   aab53:   48 8b 42 10             mov    0x10(%rdx),%rax
   aab57:   48 8b 4a 20             mov    0x20(%rdx),%rcx
   aab5b:   4c 89 bd 68 ff ff ff    mov    %r15,-0x98(%rbp)
   aab62:   4c 8d 24 c8             lea    (%rax,%rcx,8),%r12
   aab66:   8b 87 58 02 00 00       mov    0x258(%rdi),%eax
   aab6c:   44 8d 2c 85 00 00 00    lea    0x0(,%rax,4),%r13d
   aab73:   00
   aab74:   48 8d 45 90             lea    -0x70(%rbp),%rax
   aab78:   48 89 c7                mov    %rax,%rdi
   aab7b:   48 89 85 60 ff ff ff    mov    %rax,-0xa0(%rbp)
   aab82:   e8 a9 ae ff ff          callq  a5a30 <_ZNSt6vectorIN6Nektar3AVX7VecDataIdLi4EEEN5boost9alignment17aligned_allocatorIS3_Lm32EEEEC2EmRKS7_.constprop.404>
   aab87:   48 63 b3 58 02 00 00    movslq 0x258(%rbx),%rsi
   aab8e:   48 8d 45 b0             lea    -0x50(%rbp),%rax
   aab92:   48 89 c7                mov    %rax,%rdi
   aab95:   48 89 85 58 ff ff ff    mov    %rax,-0xa8(%rbp)
   aab9c:   e8 8f ae ff ff          callq  a5a30 <_ZNSt6vectorIN6Nektar3AVX7VecDataIdLi4EEEN5boost9alignment17aligned_allocatorIS3_Lm32EEEEC2EmRKS7_.constprop.404>
   aaba1:   8b 43 30                mov    0x30(%rbx),%eax
   aaba4:   85 c0                   test   %eax,%eax
   aaba6:   0f 8e f3 02 00 00       jle    aae9f <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x38f>
   aabac:   49 63 c5                movslq %r13d,%rax
   aabaf:   c5 e1 57 db             vxorpd %xmm3,%xmm3,%xmm3
   aabb3:   48 c1 e0 03             shl    $0x3,%rax
   aabb7:   4d 89 fd                mov    %r15,%r13
   aabba:   45 31 ff                xor    %r15d,%r15d
   aabbd:   48 89 85 50 ff ff ff    mov    %rax,-0xb0(%rbp)
   aabc4:   0f 1f 40 00             nopl   0x0(%rax)
   aabc8:   4c 89 e8                mov    %r13,%rax
   aabcb:   48 2b 85 68 ff ff ff    sub    -0x98(%rbp),%rax
   aabd2:   48 8b 95 60 ff ff ff    mov    -0xa0(%rbp),%rdx
   aabd9:   48 03 83 50 01 00 00    add    0x150(%rbx),%rax
   aabe0:   be 09 00 00 00          mov    $0x9,%esi
   aabe5:   4c 89 ef                mov    %r13,%rdi
   aabe8:   c5 fd 29 9d 70 ff ff    vmovapd %ymm3,-0x90(%rbp)
   aabef:   ff
   aabf0:   49 89 c6                mov    %rax,%r14
   aabf3:   c5 f8 77                vzeroupper
   aabf6:   e8 65 c0 f7 ff          callq  26c60 <_ZN6Nektar3AVX7VecDataIdLi4EE15load_interleaveEPKdmRSt6vectorIS2_N5boost9alignment17aligned_allocatorIS2_Lm32EEEE@plt>
   aabfb:   4c 8b 53 48             mov    0x48(%rbx),%r10
   aabff:   48 8b 55 90             mov    -0x70(%rbp),%rdx
   aac03:   48 8b bb 08 01 00 00    mov    0x108(%rbx),%rdi
   aac0a:   c5 fd 28 9d 70 ff ff    vmovapd -0x90(%rbp),%ymm3
   aac11:   ff
   aac12:   48 8b 73 60             mov    0x60(%rbx),%rsi
   aac16:   48 8b 8b 20 01 00 00    mov    0x120(%rbx),%rcx
   aac1d:   c4 c1 7d 28 22          vmovapd (%r10),%ymm4
   aac22:   4c 8b 5d b0             mov    -0x50(%rbp),%r11
   aac26:   c4 c1 7d 28 52 20       vmovapd 0x20(%r10),%ymm2
   aac2c:   c5 dd 59 2a             vmulpd (%rdx),%ymm4,%ymm5
   aac30:   c5 7d 28 0f             vmovapd (%rdi),%ymm9
   aac34:   c5 ed 59 4a 20          vmulpd 0x20(%rdx),%ymm2,%ymm1
   aac39:   c5 7d 28 47 20          vmovapd 0x20(%rdi),%ymm8
   aac3e:   c4 c1 7d 28 42 40       vmovapd 0x40(%r10),%ymm0
   aac44:   c5 fd 28 7f 40          vmovapd 0x40(%rdi),%ymm7
   aac49:   c5 dd 59 72 60          vmulpd 0x60(%rdx),%ymm4,%ymm6
   aac4e:   c4 c1 55 59 2e          vmulpd (%r14),%ymm5,%ymm5
   aac53:   c4 c1 75 59 4e 20       vmulpd 0x20(%r14),%ymm1,%ymm1
   aac59:   c5 dd 59 a2 c0 00 00    vmulpd 0xc0(%rdx),%ymm4,%ymm4
   aac60:   00
   aac61:   c4 c1 4d 59 76 60       vmulpd 0x60(%r14),%ymm6,%ymm6
   aac67:   c4 c2 e5 98 e9          vfmadd132pd %ymm9,%ymm3,%ymm5
   aac6c:   c4 c2 f5 b8 e8          vfmadd231pd %ymm8,%ymm1,%ymm5
   aac71:   c5 fd 59 4a 40          vmulpd 0x40(%rdx),%ymm0,%ymm1
   aac76:   c4 c1 5d 59 a6 c0 00    vmulpd 0xc0(%r14),%ymm4,%ymm4
   aac7d:   00 00
   aac7f:   c4 c2 e5 98 f1          vfmadd132pd %ymm9,%ymm3,%ymm6
   aac84:   c4 c1 75 59 4e 40       vmulpd 0x40(%r14),%ymm1,%ymm1
   aac8a:   c4 c2 e5 98 e1          vfmadd132pd %ymm9,%ymm3,%ymm4
   aac8f:   c4 e2 d5 98 cf          vfmadd132pd %ymm7,%ymm5,%ymm1
   aac94:   c5 ed 59 aa 80 00 00    vmulpd 0x80(%rdx),%ymm2,%ymm5
   aac9b:   00
   aac9c:   c5 ed 59 92 e0 00 00    vmulpd 0xe0(%rdx),%ymm2,%ymm2
   aaca3:   00
   aaca4:   c4 c1 55 59 ae 80 00    vmulpd 0x80(%r14),%ymm5,%ymm5
   aacab:   00 00
   aacad:   c4 c1 6d 59 96 e0 00    vmulpd 0xe0(%r14),%ymm2,%ymm2
   aacb4:   00 00
   aacb6:   c4 c2 d5 b8 f0          vfmadd231pd %ymm8,%ymm5,%ymm6
   aacbb:   c5 fd 59 aa a0 00 00    vmulpd 0xa0(%rdx),%ymm0,%ymm5
   aacc2:   00
   aacc3:   c5 fd 59 82 00 01 00    vmulpd 0x100(%rdx),%ymm0,%ymm0
   aacca:   00
   aaccb:   c4 c2 dd 98 d0          vfmadd132pd %ymm8,%ymm4,%ymm2
   aacd0:   c5 f5 59 26             vmulpd (%rsi),%ymm1,%ymm4
   aacd4:   c4 c1 55 59 ae a0 00    vmulpd 0xa0(%r14),%ymm5,%ymm5
   aacdb:   00 00
   aacdd:   c4 c1 7d 59 86 00 01    vmulpd 0x100(%r14),%ymm0,%ymm0
   aace4:   00 00
   aace6:   c4 e2 e5 98 21          vfmadd132pd (%rcx),%ymm3,%ymm4
   aaceb:   c4 e2 cd 98 ef          vfmadd132pd %ymm7,%ymm6,%ymm5
   aacf0:   c4 e2 ed 98 c7          vfmadd132pd %ymm7,%ymm2,%ymm0
   aacf5:   c5 d5 59 56 20          vmulpd 0x20(%rsi),%ymm5,%ymm2
   aacfa:   c4 e2 ed b8 61 20       vfmadd231pd 0x20(%rcx),%ymm2,%ymm4
   aad00:   c5 fd 59 56 40          vmulpd 0x40(%rsi),%ymm0,%ymm2
   aad05:   c4 e2 dd 98 51 40       vfmadd132pd 0x40(%rcx),%ymm4,%ymm2
   aad0b:   c4 c1 7d 29 13          vmovapd %ymm2,(%r11)
   aad10:   c5 f5 59 4e 60          vmulpd 0x60(%rsi),%ymm1,%ymm1
   aad15:   c5 d5 59 ae 80 00 00    vmulpd 0x80(%rsi),%ymm5,%ymm5
   aad1c:   00
   aad1d:   c5 fd 59 86 a0 00 00    vmulpd 0xa0(%rsi),%ymm0,%ymm0
   aad24:   00
   aad25:   c4 e2 e5 98 09          vfmadd132pd (%rcx),%ymm3,%ymm1
   aad2a:   c4 e2 f5 98 69 20       vfmadd132pd 0x20(%rcx),%ymm1,%ymm5
   aad30:   c4 e2 d5 98 41 40       vfmadd132pd 0x40(%rcx),%ymm5,%ymm0
   aad36:   c4 c1 7d 29 43 40       vmovapd %ymm0,0x40(%r11)
   aad3c:   c4 c1 7d 28 6a 60       vmovapd 0x60(%r10),%ymm5
   aad42:   c4 c1 7d 28 a2 80 00    vmovapd 0x80(%r10),%ymm4
   aad49:   00 00
   aad4b:   c5 d5 59 0a             vmulpd (%rdx),%ymm5,%ymm1
   aad4f:   c5 7d 28 0f             vmovapd (%rdi),%ymm9
   aad53:   c5 dd 59 52 20          vmulpd 0x20(%rdx),%ymm4,%ymm2
   aad58:   c5 7d 28 47 20          vmovapd 0x20(%rdi),%ymm8
   aad5d:   c4 c1 7d 28 82 a0 00    vmovapd 0xa0(%r10),%ymm0
   aad64:   00 00
   aad66:   c5 fd 28 7f 40          vmovapd 0x40(%rdi),%ymm7
   aad6b:   c5 d5 59 72 60          vmulpd 0x60(%rdx),%ymm5,%ymm6
   aad70:   c4 c1 75 59 0e          vmulpd (%r14),%ymm1,%ymm1
   aad75:   c4 c1 6d 59 56 20       vmulpd 0x20(%r14),%ymm2,%ymm2
   aad7b:   c5 d5 59 aa c0 00 00    vmulpd 0xc0(%rdx),%ymm5,%ymm5
   aad82:   00
   aad83:   c4 c1 4d 59 76 60       vmulpd 0x60(%r14),%ymm6,%ymm6
   aad89:   c4 c2 e5 98 c9          vfmadd132pd %ymm9,%ymm3,%ymm1
   aad8e:   c4 c2 ed b8 c8          vfmadd231pd %ymm8,%ymm2,%ymm1
   aad93:   c5 fd 59 52 40          vmulpd 0x40(%rdx),%ymm0,%ymm2
   aad98:   c4 c2 e5 98 f1          vfmadd132pd %ymm9,%ymm3,%ymm6
   aad9d:   c4 c1 6d 59 56 40       vmulpd 0x40(%r14),%ymm2,%ymm2
   aada3:   c4 e2 f5 98 d7          vfmadd132pd %ymm7,%ymm1,%ymm2
   aada8:   c5 dd 59 8a 80 00 00    vmulpd 0x80(%rdx),%ymm4,%ymm1
   aadaf:   00
   aadb0:   c4 c1 75 59 8e 80 00    vmulpd 0x80(%r14),%ymm1,%ymm1
   aadb7:   00 00
   aadb9:   c4 c2 f5 b8 f0          vfmadd231pd %ymm8,%ymm1,%ymm6
   aadbe:   c5 fd 59 8a a0 00 00    vmulpd 0xa0(%rdx),%ymm0,%ymm1
   aadc5:   00
   aadc6:   c4 c1 75 59 8e a0 00    vmulpd 0xa0(%r14),%ymm1,%ymm1
   aadcd:   00 00
   aadcf:   c4 c1 55 59 ae c0 00    vmulpd 0xc0(%r14),%ymm5,%ymm5
   aadd6:   00 00
   aadd8:   48 8b bd 58 ff ff ff    mov    -0xa8(%rbp),%rdi
   aaddf:   c5 dd 59 a2 e0 00 00    vmulpd 0xe0(%rdx),%ymm4,%ymm4
   aade6:   00
   aade7:   c5 fd 59 82 00 01 00    vmulpd 0x100(%rdx),%ymm0,%ymm0
   aadee:   00
   aadef:   4c 89 e2                mov    %r12,%rdx
   aadf2:   c4 e2 cd 98 cf          vfmadd132pd %ymm7,%ymm6,%ymm1
   aadf7:   c4 c2 e5 98 e9          vfmadd132pd %ymm9,%ymm3,%ymm5
   aadfc:   c4 c1 5d 59 a6 e0 00    vmulpd 0xe0(%r14),%ymm4,%ymm4
   aae03:   00 00
   aae05:   c4 c1 7d 59 86 00 01    vmulpd 0x100(%r14),%ymm0,%ymm0
   aae0c:   00 00
   aae0e:   c4 c2 d5 98 e0          vfmadd132pd %ymm8,%ymm5,%ymm4
   aae13:   c5 ed 59 2e             vmulpd (%rsi),%ymm2,%ymm5
   aae17:   c4 e2 dd 98 c7          vfmadd132pd %ymm7,%ymm4,%ymm0
   aae1c:   c5 f5 59 66 20          vmulpd 0x20(%rsi),%ymm1,%ymm4
   aae21:   c4 e2 e5 98 29          vfmadd132pd (%rcx),%ymm3,%ymm5
   aae26:   c4 e2 dd b8 69 20       vfmadd231pd 0x20(%rcx),%ymm4,%ymm5
   aae2c:   c5 fd 59 66 40          vmulpd 0x40(%rsi),%ymm0,%ymm4
   aae31:   c4 e2 d5 98 61 40       vfmadd132pd 0x40(%rcx),%ymm5,%ymm4
   aae37:   c4 c1 7d 29 63 20       vmovapd %ymm4,0x20(%r11)
   aae3d:   c5 ed 59 56 60          vmulpd 0x60(%rsi),%ymm2,%ymm2
   aae42:   c5 f5 59 8e 80 00 00    vmulpd 0x80(%rsi),%ymm1,%ymm1
   aae49:   00
   aae4a:   c5 fd 59 86 a0 00 00    vmulpd 0xa0(%rsi),%ymm0,%ymm0
   aae51:   00
   aae52:   48 63 b3 58 02 00 00    movslq 0x258(%rbx),%rsi
   aae59:   c4 e2 e5 98 11          vfmadd132pd (%rcx),%ymm3,%ymm2
   aae5e:   c4 e2 ed 98 49 20       vfmadd132pd 0x20(%rcx),%ymm2,%ymm1
   aae64:   c4 e2 f5 98 41 40       vfmadd132pd 0x40(%rcx),%ymm1,%ymm0
   aae6a:   c4 c1 7d 29 43 60       vmovapd %ymm0,0x60(%r11)
   aae70:   c5 f8 77                vzeroupper
   aae73:   e8 f8 c3 f7 ff          callq  27270 <_ZN6Nektar3AVX7VecDataIdLi4EE18deinterleave_storeERKSt6vectorIS2_N5boost9alignment17aligned_allocatorIS2_Lm32EEEEmPd@plt>
   aae78:   49 81 c5 20 01 00 00    add    $0x120,%r13
   aae7f:   4c 03 a5 50 ff ff ff    add    -0xb0(%rbp),%r12
   aae86:   41 83 c7 01             add    $0x1,%r15d
   aae8a:   44 39 7b 30             cmp    %r15d,0x30(%rbx)
   aae8e:   c5 fd 28 9d 70 ff ff    vmovapd -0x90(%rbp),%ymm3
   aae95:   ff
   aae96:   0f 8f 2c fd ff ff       jg     aabc8 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0xb8>
   aae9c:   c5 f8 77                vzeroupper
   aae9f:   48 8b 7d b0             mov    -0x50(%rbp),%rdi
   aaea3:   48 85 ff                test   %rdi,%rdi
   aaea6:   74 05                   je     aaead <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x39d>
   aaea8:   e8 43 bc f7 ff          callq  26af0 <free@plt>
   aaead:   48 8b 7d 90             mov    -0x70(%rbp),%rdi
   aaeb1:   48 85 ff                test   %rdi,%rdi
   aaeb4:   74 05                   je     aaebb <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x3ab>
   aaeb6:   e8 35 bc f7 ff          callq  26af0 <free@plt>
   aaebb:   48 8b 45 c8             mov    -0x38(%rbp),%rax
   aaebf:   64 48 33 04 25 28 00    xor    %fs:0x28,%rax
   aaec6:   00 00
   aaec8:   75 31                   jne    aaefb <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x3eb>
   aaeca:   48 83 ec 80             sub    $0xffffffffffffff80,%rsp
   aaece:   5b                      pop    %rbx
   aaecf:   41 5a                   pop    %r10
   aaed1:   41 5c                   pop    %r12
   aaed3:   41 5d                   pop    %r13
   aaed5:   41 5e                   pop    %r14
   aaed7:   41 5f                   pop    %r15
   aaed9:   5d                      pop    %rbp
   aaeda:   49 8d 62 f8             lea    -0x8(%r10),%rsp
   aaede:   c3                      retq
   aaedf:   48 8b 7d 90             mov    -0x70(%rbp),%rdi
   aaee3:   48 89 c3                mov    %rax,%rbx
   aaee6:   48 85 ff                test   %rdi,%rdi
   aaee9:   74 15                   je     aaf00 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x3f0>
   aaeeb:   c5 f8 77                vzeroupper
   aaeee:   e8 fd bb f7 ff          callq  26af0 <free@plt>
   aaef3:   48 89 df                mov    %rbx,%rdi
   aaef6:   e8 45 c3 f7 ff          callq  27240 <_Unwind_Resume@plt>
   aaefb:   e8 70 ba f7 ff          callq  26970 <__stack_chk_fail@plt>
   aaf00:   c5 f8 77                vzeroupper
   aaf03:   eb ee                   jmp    aaef3 <_ZN6Nektar3AVX15AVXIProductQuadILi4ELb1EE19AVXIProductQuadImplILi2ELi2ELi3ELi3EEEvRKNS_5ArrayINS_4OneDEKdEERNS4_IS5_dEE+0x3e3>
   aaf05:   66 2e 0f 1f 84 00 00    nopw   %cs:0x0(%rax,%rax,1)
   aaf0c:   00 00 00
   aaf0f:   90                      nop