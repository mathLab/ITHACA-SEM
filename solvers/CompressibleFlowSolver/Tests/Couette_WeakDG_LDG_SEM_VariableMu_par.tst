<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, Variable Viscosity, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Couette_WeakDG_LDG_SEM_VariableMu_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_VariableMu_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-7">0.0805614</value>
            <value variable="rhou" tolerance="1e-4">51.8867</value>
            <value variable="rhov" tolerance="1e-6">0.222213</value>
            <value variable="E" tolerance="1e-2">4415.64</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0716895</value>
            <value variable="rhou" tolerance="1e-4">55.3017</value>
            <value variable="rhov" tolerance="1e-6  ">0.424398</value>
            <value variable="E" tolerance="1e-2">4227.23</value>
        </metric>
    </metrics>
</test>


