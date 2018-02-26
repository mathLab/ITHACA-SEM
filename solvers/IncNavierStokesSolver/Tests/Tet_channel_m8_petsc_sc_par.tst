<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Tet elements, PETSc sc, par(3), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch -I GlobalSysSoln=PETScStaticCond Tet_channel_m8_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Tet_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">8.62254e-15</value>
            <value variable="v" tolerance="1e-12">7.34883e-15</value>
            <value variable="w" tolerance="1e-12">3.03011e-14</value>
            <value variable="p" tolerance="1e-8">4.64664e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.70349e-14</value>
            <value variable="v" tolerance="1e-12">2.65691e-14</value>
            <value variable="w" tolerance="1e-12">1.05027e-13</value>
            <value variable="p" tolerance="1e-8">1.40776e-12</value>
        </metric>
    </metrics>
</test>
