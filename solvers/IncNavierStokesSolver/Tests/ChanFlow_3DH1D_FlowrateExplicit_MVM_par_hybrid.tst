<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=4, 12 Fourier modes (MVM) with flowrate control and forcing direction explicitly defined</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--npz 2 ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">2.59003e-08</value>
            <value variable="v" tolerance="1e-8">6.79495e-11</value>
            <value variable="w" tolerance="1e-8">0</value>
            <value variable="p" tolerance="1e-6">1.95738e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">4.67694e-08</value>
            <value variable="v" tolerance="1e-8">2.2142e-10</value>
            <value variable="w" tolerance="1e-8">5.01397e-19</value>
            <value variable="p" tolerance="1e-6">3.62234e-07</value>
        </metric>
    </metrics>
</test>
