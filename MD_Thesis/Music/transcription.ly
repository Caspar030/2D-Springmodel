\version "2.24.0"

\header {
  arranger = "Tunescribers"
  tagline  = ##f
}

chordsA = \chordmode {
  \time 5/4
  ees1*5/4 |
  \time 4/4
  g4:m f4 c2:7 |
  f1 |
  bes2 aes2:maj9 |
  \time 5/4
  bes1:7 aes4:maj7 |
  bes1:7 aes4:maj7 |
  bes1*5/4:7
}

melody = \relative c' {
  \time 5/4
  % --- m1 : Eb ---
  g8 b8 a8 c8. c16 g8 b8 a8 c8. c16 |
  \time 4/4
  % --- m2 : Gm F C7 ---
  g8 a8 bes8 b8 c4 c4 |
  % --- m3 : F ---
  b'2~ b4 c8 d8 |
  % --- m4 : Bb  Abmaj9 ---
  bes4 d4 <aes c' ees' g'>2 |
  \time 5/4
  % --- m5 : Bb7  Abmaj7 ---
  bes8 d8 f8 aes8. aes16 aes8 c8 ees8 g8. g16 |
  % --- m6 : Bb7  Abmaj7 (repeat of m5) ---
  bes8 d8 f8 aes8. aes16 aes8 c8 ees8 g8. g16 |
  % --- m7 : Bb7 (final) ---
  bes8 d8 f8 aes8 g8 f8 ees2
}

bassLine = \relative c {
  \time 5/4
  % --- m1 ---
  ees8 bes8 bes8 ees8 bes8 ees8 bes8 bes8 ees8 bes8 |
  \time 4/4
  % --- m2 ---
  g4 f4 c4 c4 |
  % --- m3 ---
  f2 bes4 f4 |
  % --- m4 ---
  bes4 f4 aes2 |
  \time 5/4
  % --- m5 ---
  bes8 f8 bes8 f8 bes8 aes8 ees8 aes8 ees8 aes8 |
  % --- m6 (repeat of m5) ---
  bes8 f8 bes8 f8 bes8 aes8 ees8 aes8 ees8 aes8 |
  % --- m7 ---
  bes8 f8~ f2~ f4. bes8
}

\score {
  <<
    \new ChordNames \chordsA
    \new PianoStaff <<
      \new Staff = "right" \with { instrumentName = "Piano" }
      {
        \key ees \major
        \tempo 4 = 92
        \mark "A"
        \clef treble
        \melody
      }
      \new Staff = "left"
      {
        \key ees \major
        \clef bass
        \bassLine
      }
    >>
  >>
  \layout { }
  \midi { }
}
