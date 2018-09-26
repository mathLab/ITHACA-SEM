<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=5, 8 Fourier modes (MVM) with flowrate control</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_Flowrate_MVM.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_Flowrate_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">6.57092e-15</value>
            <value variable="p" tolerance="1e-8">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.93546e-18</value>
            <value variable="v" tolerance="1e-12">2.93958e-18</value>
            <value variable="w" tolerance="1e-12">2.37588e-14</value>
            <value variable="p" tolerance="1e-8">2.74618e-17</value>
        </metric>
    </metrics>
</test>
