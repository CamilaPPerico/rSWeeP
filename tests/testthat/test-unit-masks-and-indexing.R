# test-unit-masks-and-indexing.R: Etapas 1 e 2 do Artigo (Spaced Words e Indexação Vigesimal/Quaternária)

test_that("convertMask e Defmask convertem máscaras em formato binário expandido", {
  # Máscara canônica do artigo para proteínas: c(2, 1, 2) -> 11011
  expect_equal(rSWeeP:::convertMask(c(2, 1, 2)), c(1, 1, 0, 1, 1))
  
  # Máscara para nucleotídeos: c(5, 5, 5) -> 5 uns, 5 zeros, 5 uns
  expect_equal(rSWeeP:::convertMask(c(5, 5, 5)), c(rep(1, 5), rep(0, 5), rep(1, 5)))

  # Defmask com máscara default
  expect_equal(rSWeeP:::Defmask(NULL, "AA"), c(1, 1, 0, 1, 1))
  expect_equal(rSWeeP:::Defmask(NULL, "NT"), c(rep(1, 5), rep(0, 5), rep(1, 5)))

  # Defmask com máscara já binária
  expect_equal(rSWeeP:::Defmask(c(1, 0, 1), "AA"), c(1, 0, 1))
  
  # Defmask: se max(mask) <= 1, mantém inalterado
  expect_equal(rSWeeP:::Defmask(c(1, 1, 1), "AA"), c(1, 1, 1))
})

test_that("seq2num mapeia alfabeto canônico e caracteres especiais para inteiros", {
  # Alfabeto de 20 aminoácidos canônicos -> 1:20
  aa_canon <- "ARNDCQEGHILKMFPSTWYV"
  expect_equal(rSWeeP:::seq2num(aa_canon, "AA"), 1:20)

  # Caracteres ambíguos, raros e especiais de aminoácidos
  # B=21, Z=22, X=23, *=24, -=25, U=26, O=27, J=28, ?=0
  aa_special <- "BZX*-UOJ?"
  expect_equal(rSWeeP:::seq2num(aa_special, "AA"), c(21, 22, 23, 24, 25, 26, 27, 28, 0))

  # Tolerância a caixa baixa
  expect_equal(rSWeeP:::seq2num("arnd", "AA"), 1:4)

  # Alfabeto canônico e ambiguidades de nucleotídeos
  # A=1, C=2, G=3, T=4, U=4, R=5, Y=6, K=7, M=8, S=9, W=10, B=11, D=12, H=13, V=14, N=15, -=16, *=0
  nt_test <- "ACGTU"
  expect_equal(rSWeeP:::seq2num(nt_test, "NT"), c(1, 2, 3, 4, 4))
  expect_equal(rSWeeP:::seq2num("N*", "NT"), c(15, 0))
})

test_that("returnPOT calcula índices de coocorrência vigesimal e quaternária bijetiva", {
  # Proteína: máscara c(1, 0, 1) sobre "ACD" -> k-mers "A.D"
  # A=1 (val-1 = 0), C=2, D=4 (val-1 = 3)
  # Indexação f(S_k) = 1 + sum((V(S_q) - 1) * 20^(q-1))
  # Para "A.D": pos1=A (0), pos2=D (3) -> 1 + 0*1 + 3*20 = 61
  seq_num <- rSWeeP:::seq2num("ACD", "AA")
  expect_equal(rSWeeP:::returnPOT(seq_num, c(1, 0, 1), "AA"), 61)

  # Sequência homogênea "AAAA" com máscara 11 (k=2) -> 1 + 0 + 0 = 1
  seq_aaaa <- rSWeeP:::seq2num("AAAA", "AA")
  expect_equal(rSWeeP:::returnPOT(seq_aaaa, c(1, 1), "AA"), c(1, 1, 1))

  # Sequência com resíduo máximo "VV" com máscara 11 (k=2, V=20, val-1=19) -> 1 + 19 + 19*20 = 400
  seq_vv <- rSWeeP:::seq2num("VV", "AA")
  expect_equal(rSWeeP:::returnPOT(seq_vv, c(1, 1), "AA"), 400)

  # Delimitador '*' (val=24 > 19) e '?' (val=0 < 1) devem ser ignorados (pot <= 0 filtrado)
  seq_with_spacer <- rSWeeP:::seq2num("AC*DE", "AA")
  # Janelas que cruzam '*' em posições ativas recebem pot = -1 e não são retornadas.
  # "AC*" e "*DE" são invalidadas; apenas "C*D" sobrevive (C=5, D=4 -> 1 + 4*1 + 3*20 = 65)
  pots <- rSWeeP:::returnPOT(seq_with_spacer, c(1, 0, 1), "AA")
  expect_equal(pots, 65)
})
