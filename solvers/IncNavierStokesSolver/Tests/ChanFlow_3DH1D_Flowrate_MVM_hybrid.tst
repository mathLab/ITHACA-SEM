<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=5, 8 Fourier modes (MVM) with flowrate control</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--npz 2 ChanFlow_3DH1D_Flowrate_MVM.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_Flowrate_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.54851e-15</value>
            <value variable="v" tolerance="1e-12">1.92039e-15</value>
            <value variable="w" tolerance="1e-12">8.94953e-07</value>
            <value variable="p" tolerance="1e-8">3.23051e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">9.74693e-15</value>
            <value variable="v" tolerance="1e-12">7.87791e-15</value>
            <value variable="w" tolerance="1e-12">2.77108e-06</value>
            <value variable="p" tolerance="1e-8">1.64691e-14</value>
        </metric>
    </metrics>
</test>
