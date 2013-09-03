import EclipseKeys._

lazy val Jandom = project in file("core")

lazy val JandomExtended = project in file("extended") dependsOn Jandom % "compile->compile;test->test"

lazy val root = project in file(".") aggregate (Jandom, JandomExtended) 

version in ThisBuild := "0.1.2-SNAPSHOT"

scalaVersion in ThisBuild := "2.10.2"

executionEnvironment in ThisBuild := Some(EclipseExecutionEnvironment.JavaSE17)

fork in ThisBuild := true

// This delegates the root run task to the run task in JandomCore

run <<= run in ("Jandom", Compile)

