bool IsMIDTrackletSelected(TVector3 posHitLayer1, TVector3 posHitLayer2, ) {

  TVector3 primVtx(0,0,0);

  return IsMIDTrackletSelected(posHitLayer1,posHitLayer2,primVtx);

}

//====================================================================================================================================================

bool IsMIDTrackletSelected(TVector3 posHitLayer1, TVector3 posHitLayer2, TVector3 primVtx) {

  // Draft, to be finalized
  
  double deltaPhi = posHitLayer1.Phi() - posHitLayer2.Phi();
  double deltaEta = posHitLayer1.Eta() - posHitLayer2.Eta();
  
  if (abs(deltaPhi) > 0.06) return kFALSE;
  if (abs(deltaEta) > 0.06) return kFALSE;

  return kTRUE;

}

//====================================================================================================================================================

bool IsMIDTrackletSelected(TVector3 posHitLayer1, TVector3 posHitLayer2, TVector3 primVtx, TVector3 trackITS) {

  // Draft, to be finalized
  
  double deltaPhi = posHitLayer1.Phi() - posHitLayer2.Phi();
  double deltaEta = posHitLayer1.Eta() - posHitLayer2.Eta();
  
  if (abs(deltaPhi) > 0.06) return kFALSE;
  if (abs(deltaEta) > 0.06) return kFALSE;

  return kTRUE;

}

//====================================================================================================================================================

