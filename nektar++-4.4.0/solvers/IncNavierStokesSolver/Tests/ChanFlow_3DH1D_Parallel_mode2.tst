<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D homogeneous 1D Channel Flow, HOM parallelisation (2 proc)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--npz 2 ChanFlow_3DH1D_Parallel_mode2.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_Parallel_mode2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">3.57865e-14</value>
            <value variable="v" tolerance="1e-6">2.17691e-13</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">5.03999e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">2.83162e-13</value>
            <value variable="v" tolerance="1e-6">3.42167e-13</value>
            <value variable="w" tolerance="1e-6">1.22076e-18</value>
            <value variable="p" tolerance="1e-6">1.39357e-10</value>
        </metric>
    </metrics>
</test>


