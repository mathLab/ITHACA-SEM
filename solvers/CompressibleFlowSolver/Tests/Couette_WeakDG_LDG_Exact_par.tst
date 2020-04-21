<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with periodic BCs, exact solution, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>-I parts=6,7:2,4,5:0,1,3  Couette_WeakDG_LDG_Exact_par.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_Exact_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">5.26e-11</value>
            <value variable="rhou" tolerance="1e-12">1.06814e-09</value>
            <value variable="rhov" tolerance="1e-12">1.13795e-09</value>
            <value variable="E" tolerance="1e-8">1.06807e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.72468e-10</value>
            <value variable="rhou" tolerance="1e-12">1.30154e-09</value>
            <value variable="rhov" tolerance="1e-12">1.03028e-09</value>
            <value variable="E" tolerance="1e-8">8.58388e-07</value>
        </metric>
    </metrics>
</test>
